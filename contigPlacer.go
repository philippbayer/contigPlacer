/*
This is the source-code for contigPlacer, Philipp Bayer 2013, MIT-license.

contigPlacer parses a gff3-file of already placed contigs and compares these to a map and a dat-file of SNPs from Flapjack, placing all SNPs on their respective contigs in the process, and combining all SNPs into one metaSNP. It then parses a gff3-file of unplaced contigs (and a map and a dat-file of SNPs) to create a list of unplaced contigs with their respective metaSNPs.

Then, all unplaced contigs are compared against all placed contigs using a penalized Hamming distance, and accordingly placed. A new gff3-file is produced showing where all unplaced contigs and placed contigs are in relation to each other on which chromosome, one file for each chromosome.

Usage:
Usage of ./contigPlacer:
-chrom="": The path to the file detailing assembled chromosomes.
-contig="": The path to the file detailing contigs.
-v=false: Prints the current version.
-s: Writes the best scores for each unplaced contig into "all_scores.txt"

Format of an input chrom-file (space-delimited) is generally "dat-file map-file gff3-file", absolute or relative paths allowed.
GBS_tapidor.chrA01_fixed_cleaned_imputed_errorindividuals_removed.dat GBS_tapidor.chrA01_fixed_cleaned.map tapidor_chrA01_pseudo_contig.gff3
GBS_tapidor.chrA02_fixed_cleaned_imputed_errorindividuals_removed.dat GBS_tapidor.chrA02_fixed_cleaned.map tapidor_chrA02_pseudo_contig.gff3
GBS_tapidor.chrA03_fixed_cleaned_imputed_errorindividuals_removed.dat GBS_tapidor.chrA03_fixed_cleaned.map tapidor_chrA03_pseudo_contig.gff3
GBS_tapidor.chrA04_fixed_cleaned_imputed_errorindividuals_removed.dat GBS_tapidor.chrA04_fixed_cleaned.map tapidor_chrA04_pseudo_contig.gff3

Format of an input contig file:
GBS_tapidor.chrA01.to.be.placed_fixed_cleaned_imputed_errorindividuals_removed.dat GBS_tapidor.chrA01.to.be.placed_fixed_cleaned.map tapidor_chrA01_to_be_placed_pseudo_contig.gff3
GBS_tapidor.chrA02.to.be.placed_fixed_cleaned_imputed_errorindividuals_removed.dat GBS_tapidor.chrA02.to.be.placed_fixed_cleaned.map tapidor_chrA02_to_be_placed_pseudo_contig.gff3
GBS_tapidor.chrA03.to.be.placed_fixed_cleaned_imputed_errorindividuals_removed.dat GBS_tapidor.chrA03.to.be.placed_fixed_cleaned.map tapidor_chrA03_to_be_placed_pseudo_contig.gff3

TODO:
      1. Implement concurrency - huge speed gains, but complicated
      2. Implement ability to send e-mails.

Tested with:
Ubuntu 12.04, Fedora 18, 22, Go 1.1 to 1.2.1, 1.5.1
*/

package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
	"sort"
	"strconv"
	"strings"
)

const (
	SNPNameDelimiter string  = "|||"                  // A hackish delimiter so we can store all SNPs into one map without (hopefully!) collisions
	SNPMinimumRange  int     = 1                      // SNPs closer than this to other already stored SNPs are not included
	CurrentVersion   float64 = 2.1                    // Because Torsten Seemann said so.
	OutputSuffix     string  = "_contigs_placed.gff3" // The suffix of the output files written by this program
	ScoresFile       string  = "all_scores.txt"       // The name of the file with all scores
	UnplaceableFile  string  = "list_of_unplaceable_contigs.txt"
)

var SkippableThings = map[string]string{"": "", " ": "", "E": "", "e": ""} // These are alleles we ignore when comparing two SNPs
var ComparisonMethods = map[string]string{"hamming": "", "r_squared": "", "percentage": ""}

// It's not possible to declare maps as constants. :(

//##############################################################################################################################

type Contig struct {
	name, chromosome, orientation string
	startPos, endPos              int
	nextContig                    *Contig
	metaSNP, lastMeta, listOfSNPs []string              // ["N", "T", "N", "N", "", ...]
	dictOfPartners                map[float64][]*Contig // Only used in newly placed contigs - key: score of similarity, value: slice of unplaced contig
	normalisedRecombinations      float64               // Stores how often we see double or more recombinations on this contig, divided by number of individuals
}

// A ContigMap is a simple linked list based on a gff3-file
type ContigMap struct {
	firstContig *Contig
}

//##############################################################################################################################

// Some basic check functions to see whether an element is contained in a map
func isIntInMap(element int, a_map map[int]string) (contained bool) {
	_, ok := a_map[element]
	if ok {
		return true
	}
	return false
}
func isStrInMap(element string, a_map map[string]string) (contained bool) {
	_, ok := a_map[element]
	if ok {
		return true
	}
	return false
}

//##############################################################################################################################

// Gets the name of the scaffold from a gff3-entry, example:
// ID=scaffold128127|size365;Name=scaffold128127|size365 becomes scaffold128127|size365
func getName(toCheck []string) (name string) {
	lastElement := toCheck[len(toCheck)-1]
	list := strings.Split(lastElement, ";")
	for _, element := range list {
		if strings.Count(element, "Name") > 0 {
			name = strings.Replace(element, "Name=", "", -1)
			return name
		}
	}
	log.Fatal("Couldn't find any contig names in the gff3-file. Check whether you have a 'Name=' field in your gff3-files.")
	return
}

//##############################################################################################################################

// Compares the alleles of two SNPs and returns a score based on how similar they are
func compareSNPs(first_snp []string, second_snp []string) (score float64) {
	var len_counter float64
	for index, first_allele := range first_snp {
		partner_allele := second_snp[index]
		if isStrInMap(first_allele, SkippableThings) || isStrInMap(partner_allele, SkippableThings) {
			// We don't count empty or E
			continue
		}
		if first_allele == partner_allele {
			score += 1
		}
		len_counter++
	}
	score = score / len_counter
	if len_counter == 0 {
		// Happens when one of both partners is entirely empty
		score = 0
	}

	return score // Score has to be between 0 and 1, where 1 is identical
}

//##############################################################################################################################

// Compares two SNPs using a penalized Hamming distance
// Scores:
// SNP 1 has A, SNP 2 has B, add 1
// SNP 1 has "", SNP 2 has B, and vice-versa, add 0.5 (it's after all a 50/50 chance)
// SNP 1 has "", SNP 2 has "", add 0.75, since 0.5 + 0.5**2
func compareSNPsHamming(first_snp []string, second_snp []string) (score float64) {
	if len(first_snp) != len(second_snp) {
		log.Fatal("FATAL ERROR: Two SNPs are not of identical length")
	}
	for index, first_allele := range first_snp {
		partner_allele := second_snp[index]
		if isStrInMap(first_allele, SkippableThings) && isStrInMap(partner_allele, SkippableThings) {
			// Are both empty?
			score += 0.75
		} else if isStrInMap(first_allele, SkippableThings) || isStrInMap(partner_allele, SkippableThings) {
			// Is only one empty?
			score += 0.5
		} else if first_allele != partner_allele {
			// Both are different
			score += 1
		}
	}
	return score
}

//##############################################################################################################################

// Compares alleles of two SNPs using LD, r**2
func compareSNPsLD(first_snp []string, second_snp []string) (score float64) {
	var a1, a2, b1, b2 string
	var a1b1, a1b2, a2b1, a2b2, length float64

	for index, first_allele := range first_snp {
		partner_allele := second_snp[index]
		if isStrInMap(first_allele, SkippableThings) || isStrInMap(partner_allele, SkippableThings) {
			// We don't count empty or E
			continue
		}
		if a1 == "" {
			a1 = first_allele
			b1 = first_allele
		}
		if a1 != "" && a2 == "" {
			if first_allele == a1 && partner_allele != first_allele {
				a2 = partner_allele
				b2 = partner_allele
			} else if first_allele != a1 && partner_allele != first_allele {
				a2 = first_allele
				b2 = first_allele
			}
		}
		if first_allele == a1 && partner_allele == b1 {
			a1b1 += 1
		} else if first_allele == a2 && partner_allele == b2 {
			// partner_allele is now b2
			a2b2 += 1
		} else if first_allele == a2 && partner_allele == b1 {
			a2b1 += 1
		} else if first_allele == a1 && partner_allele == b2 {
			a1b2 += 1
		}
		length += 1
	}
	if length == 0 {
		// Everything is empty.
		return 0.0
	}

	if a2 == "" || b2 == "" {
		// There is no a2, or b2 - that means all alleles are identical, that means r² == 1
		// That means we compare two monomorphic SNPs, which should never happen anyway.
		return 1
	}

	a1b1 = a1b1 / length
	a1b2 = a1b2 / length
	a2b1 = a2b1 / length
	a2b2 = a2b2 / length
	A1 := a1b1 + a1b2
	A2 := a2b1 + a2b2
	B1 := a1b1 + a2b1
	B2 := a1b2 + a2b2
	D := a1b1 - (A1 * B1)
	if D == 0 {
		return 0
	}
	r := D / (math.Sqrt(A1 * A2 * B1 * B2))
	score = math.Pow(r, 2)

	return score
}

//##############################################################################################################################

// Iterates over a slice of strings, counts how often each string appears, returns the string that appears the most
// Change in 2.1: Also counts how often a recombination (a change from A to B or B to A) occurs
func getMajorType(toCount []string) (biggestType string, recombinationCounter int) {
	var (
		biggestCount   int
		currentElement string
	)
	counts := make(map[string]int)
	for _, element := range toCount {
		if isStrInMap(element, SkippableThings) {
			continue
		}
		if element != currentElement {
			if currentElement != "" {
				recombinationCounter += 1
			}
			currentElement = element
		}

		_, ok := counts[element]
		if !ok {
			counts[element] = 1
		} else {
			counts[element] += 1
		}
	}
	for key, value := range counts {
		if value > biggestCount {
			biggestCount = value
			biggestType = key
		}
	}
	return biggestType, recombinationCounter
}

//##############################################################################################################################

// Tiny meta-function that makes the meta-SNP for a contig
// New in 2.1: also computes the number of recombinations, divided by the number of individuals
func makeMetaSNP(snps []string, HashAlleles map[string][]string) (meta []string, normalisedRecombinations float64) {
	var (
		counter           int
		allRecombinations int
	)
	length_of_alleles := len(HashAlleles[snps[0]]) // Identical to number of individuals
	meta = make([]string, length_of_alleles)
	// Iterate over each position of the SNP, and iterate over all individuals for that position
	for counter < length_of_alleles {
		alleles_on_this_position := make([]string, len(snps))
		for index, snp := range snps {
			all_alleles, ok := HashAlleles[snp]
			if !ok {
				// SNP is not in HashAlleles, was removed earlier
				continue
			}
			this_allele := strings.ToUpper(all_alleles[counter])
			if isStrInMap(this_allele, SkippableThings) {
				alleles_on_this_position[index] = ""
			} else {
				alleles_on_this_position[index] = this_allele
			}
		}
		// Now count how often each type appears, take the one that appears the most, if identical, choose randomly
		major, recombinationCounter := getMajorType(alleles_on_this_position)
		allRecombinations += recombinationCounter
		meta[counter] = major
		counter++
	}
	normalisedRecombinations = float64(allRecombinations) / float64(length_of_alleles)
	return meta, normalisedRecombinations
}

//##############################################################################################################################

// Checks gff and map files to see whether both have the same names for chromosomes, dies if different
func checkForProblems(maps []string, gffs []string) {
	var (
		file           *os.File
		err            error
		expected_chrom string
	)
	// Iterate over each pair of map- and gff-files
	// (don't need to iterate over dat-files as no chromosome-names in there)
	for index, thisMap := range maps {
		thisGff := gffs[index]

		// Get the first non-commented line in the gff3
		if file, err = os.Open(thisGff); err != nil {
			log.Fatal(err)
		}
		defer file.Close()

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			line := scanner.Text()
			if strings.Contains(line, "#") {
				continue
			}
			ll := strings.Split(line, "\t")
			// First non-comment line contains the name of the chromosome
			expected_chrom = ll[0]
			break
		}

		// Get the name from the map-file
		if file, err = os.Open(thisMap); err != nil {
			log.Fatal(err)
		}
		defer file.Close()
		scanner = bufio.NewScanner(file)
		scanner.Scan() // # FJFILE = MAP
		scanner.Scan() // Get first line
		line := strings.Split(scanner.Text(), "\t")
		map_chrom := line[1]

		// Die if not identical
		if map_chrom != expected_chrom {
			log.Fatal("FATAL ERROR: The map-file '", thisMap, "' contains the chromosome name '", map_chrom, "' but the gff-file '", thisGff, "' has '", expected_chrom, "'. Both need to be identical, please fix this.")
		}
	}
}

//##############################################################################################################################

// Helper-function - checks whether all dat-files have the same amount of lines, i.e., all SNPs have the same amount of alleles
// Iterates over all dat-files, runs "wc -l", dies if different
func checkDatFiles(list_of_dats []string) {
	length := -1
	original_dat := ""
	for _, element := range list_of_dats {
		file, err := os.Open(element)
		var out bytes.Buffer
		if err != nil {
			log.Fatal("FATAL ERROR: File ", file, " doesn't exist")
		}
		defer file.Close()
		cmd := exec.Command("wc", "-l")
		cmd.Stdin = file
		cmd.Stdout = &out
		err = cmd.Run()
		if err != nil {
			log.Fatal(err)
		}
		thisLength, err := strconv.Atoi(strings.Replace(out.String(), "\n", "", -1))
		if err != nil {
			log.Fatal(err)
		}
		if length == -1 {
			// First dat-file - nothing to compare against
			length = thisLength
			original_dat = element
			continue
		}

		if thisLength != length {
			log.Fatal("FATAL ERROR: dat-file ", element, " has a different amount of individuals (", thisLength, ") than dat-file ", original_dat, "(", length, "). Both need to be identical!")
		}
	}
}

//##############################################################################################################################

// For one snp, iterates over all possible placed SNPs and returns the best partner-SNP
func findBestPartner(SNPAlleles []string, chromMapDict map[string]*ContigMap, comparisonMethod *string) (bestPartner *Contig, bestScore float64, secondBestPartner *Contig) {
	if *comparisonMethod == "hamming" {
		bestScore = 1<<31 - 1
	} else {
		bestScore = -1.0
	}
	secondBestScore := -1.0
	for _, placedContigList := range chromMapDict {
		subContig := placedContigList.firstContig
		for subContig != nil {
			possiblePartner := subContig
			possiblePartnerAlleles := subContig.metaSNP
			if len(possiblePartnerAlleles) == 0 {
				subContig = subContig.nextContig
				continue
			}
			var score float64
			if *comparisonMethod == "percentage" {
				score = compareSNPs(SNPAlleles, possiblePartnerAlleles) // Use %
			} else if *comparisonMethod == "r_squared" {
				score = compareSNPsLD(SNPAlleles, possiblePartnerAlleles) // Use r²
			} else if *comparisonMethod == "hamming" {
				score = compareSNPsHamming(SNPAlleles, possiblePartnerAlleles) // Use Hamming
			}
			if *comparisonMethod == "hamming" {
				// Hamming score - the lower the score, the better the partners
				if score < bestScore {
					// Set the new best score
					bestPartner = possiblePartner
					bestScore = score
					continue // Don't check for second best partner
				}

				// No new best partner, check for second best partner
				if secondBestScore == -1.0 {
					// Initialize secondBestScore to first seen score
					secondBestScore = score
				} else if score > bestScore && score < secondBestScore {
					secondBestPartner = possiblePartner
					secondBestScore = score
				}
			} else {
				// All other scores, the higher the score, the better
				if score > bestScore {
					bestPartner = possiblePartner
					bestScore = score
					continue
				}
				if secondBestScore == -1.0 {
					secondBestScore = score
				} else if score < bestScore && score > secondBestScore {
					secondBestPartner = possiblePartner
					secondBestScore = score
				}
			}
			subContig = subContig.nextContig
		}
	}
	return bestPartner, bestScore, secondBestPartner
}

//##############################################################################################################################

// parses a file of gff3-names, parses each of these names
// For each name, creates a linked list of placed contigs for each gff3
func parseGffControlFile(path string, partnerMap map[string][]int, partnerChrom map[string][]string, chromMapDict map[string]*ContigMap, placedHashAlleles map[string][]string) {
	var (
		file       *os.File  // the gff3 file we parse
		chromosome string    // Name of the chromosome we're looking at
		contigMap  ContigMap // Linked list storing contigs and their order
		prevContig *Contig
		start, end int
		err        error
	)

	// Open the gff3-file
	if file, err = os.Open(path); err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	// Use bufio.Scanner because gff3-lines are not long
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		// Skip comments
		if strings.ContainsAny(line, "#") {
			continue
		}
		l := strings.Split(line, "\t")
		if start, err = strconv.Atoi(l[3]); err != nil {
			log.Fatal(err)
		}

		if end, err = strconv.Atoi(l[4]); err != nil {
			log.Fatal(err)
		}
		// gff3 is 1-inclusive, but Go and flapjack-maps are 0-inclusive, change so I don't get confused
		start -= 1
		end -= 1
		name := getName(l)

		c := new(Contig) // c is now of type *Contig
		c.startPos = start
		c.endPos = end
		c.name = name
		if chromosome == "" {
			// We have to initialize our linked list of contigs
			chromosome = l[0]
			contigMap.firstContig = c
			prevContig = c
		} else {
			prevContig.nextContig = c
			prevContig = c
		}
		c.chromosome = chromosome
		// Now get all the SNPs that are on this chromosome between start and end
		thisSNPs := make([]string, 0)
		chromNames := partnerChrom[chromosome]
		chromPos := partnerMap[chromosome]
		for index, element := range chromNames {
			position := chromPos[index]
			if position >= start && position <= end {
				thisSNPs = append(thisSNPs, element)
			} else if position > end {
				// Since the SNPs are ordered by position, we can save some time here
				break
			}
		}
		// Make the meta SNP of thisSNPs
		c.listOfSNPs = thisSNPs
		meta := []string{}
		var normalisedRecombinations float64
		if len(thisSNPs) != 0 {
			meta, normalisedRecombinations = makeMetaSNP(thisSNPs, placedHashAlleles)
		}
		c.normalisedRecombinations = normalisedRecombinations
		c.metaSNP = meta
		// Make the metaSNP of the last 5 SNPs if applicable
		lastMeta := []string{}
		if len(thisSNPs) > 10 {
			lastStart := len(thisSNPs) - len(thisSNPs)/3
			// Take at least 5 SNPs, even for smaller contigs
			if lastStart < 5 {
				lastStart = 5
			}
			lastEnd := len(thisSNPs)
			lastBit := thisSNPs[lastStart:lastEnd]
			lastMeta, _ = makeMetaSNP(lastBit, placedHashAlleles)
		}
		c.lastMeta = lastMeta
	}

	chromMapDict[chromosome] = &contigMap // Store the address of the contigMap in the big dictionary
	return
}

//##############################################################################################################################

// This function returns a slice of strings which are files to be parsed - reads entire file at once, not to use for Flapjack-files!!!
func parseControlFiles(path string) (dats []string, maps []string, gffs []string) {
	var (
		file *os.File
		err  error
	)
	// By default, the delimiter for the control-files is space
	delimiter := " "

	if file, err = os.Open(path); err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Remove any trailing white-spaces, and split by delimiter
		line := strings.Split(strings.Trim(scanner.Text(), " "), delimiter)
		if line[0] == "" || len(line) == 1 {
			// skip empty lines
			continue
		}
		if len(line) != 3 {
			log.Fatal("FATAL ERROR: Your control-file is neither tab- nor space-delimited, or doesn't include enough elements per line.")
			// Alternative: Guess other delimiters?
		}
		thisDat := line[0]
		thisMap := line[1]
		thisGff := line[2]
		dats = append(dats, thisDat)
		maps = append(maps, thisMap)
		gffs = append(gffs, thisGff)
	}
	return dats, maps, gffs
}

//##############################################################################################################################
// Tests whether bufio's Scanner can handle a file - as of 1.1.2, bufio.Scanner() dies with large files
func linesAreTooLong(f string) (tooLong bool) {
	var (
		file *os.File
		err  error
	)

	if file, err = os.Open(f); err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	scanner.Scan()
	if scanner.Text() != "# fjFile = GENOTYPE" {
		log.Fatal("FATAL ERROR: The Flapjack-dat file '" + f + "' doesn't seem to be a dat-file, header missing or wrong.\nShould be '# fjFile = GENOTYPE'.")
	}
	// With a line too long (and other problems), the next scan should trigger an error, and scanner.Err() returns an error
	scanner.Scan()
	err = scanner.Err()
	if err != nil {
		return true
	}
	return false
}

//##############################################################################################################################

// This is the simple parsing-version using bufio
func parseDatBufio(f string, hashAlleles map[string][]string, chromosome string, allSkippedSNPs map[string][]string) {
	var (
		file        *os.File
		err         error
		lineCounter int
		SNPNames    []string
	)

	indexesOfSkippedSNPs := make(map[int]string) // Since go doesn't have the "set"-type, we have to cheat slightly

	if file, err = os.Open(f); err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	skippedSNPs := allSkippedSNPs[chromosome]

	scanner := bufio.NewScanner(file)
	scanner.Scan()
	if scanner.Text() != "# fjFile = GENOTYPE" {
		log.Fatal("FATAL ERROR: The Flapjack-dat file '" + f + "' doesn't seem to be a map-file, header missing or wrong.\nShould be '# fjFile = GENOTYPE'.")
	}

	for scanner.Scan() {
		ll := strings.Split(strings.ToUpper(scanner.Text()), "\t")
		if lineCounter == 0 {
			// We have SNP names
			SNPNames = ll[1:] // In flapjack-dat-headers, the first element is empty

			// Rename SNPNames form M1 to tapidor.c09|||M1
			for i, element := range SNPNames {
				SNPNames[i] = chromosome + SNPNameDelimiter + element
			}

			// We have to get the indexes for the SNPs in SNPNames so we can skip these in later lines
			for index, element := range SNPNames {
				// Is this element in our list of snps we should skip? If yes, don't include in positions array
				include := true
				for _, toCheck := range skippedSNPs {
					if toCheck == element {
						indexesOfSkippedSNPs[index] = "y"
						include = false
						break
					}
				}
				// Initialize a new SNP in the dictionary
				if include {
					// If the SNP already exists in the dictionary, the dat-file doesn't have unique enough SNP names
					_, ok := hashAlleles[element]
					if ok {
						log.Fatal("ERROR: The SNP ", element, " exists more than once in your dat-file! Aborting.")
					}
					hashAlleles[element] = make([]string, 0)
				}
			}
		} else {
			// We have normal alleles
			ll = ll[1:]

			// In some files, the lines are one longer than the header, with an empty element - stop.
			if len(ll) != len(SNPNames) {
				log.Fatal("ERROR: The length of the line of SNP-names (line 2) is different than a line of alleles, the length of SNP names is ", len(SNPNames), ", length of alleles in line number ", lineCounter, " are ", len(ll), "\nFile is ", f)
			}

			theseSkippedSNPs := make([]string, 0)
			for index, element := range ll {
				if isIntInMap(index, indexesOfSkippedSNPs) {
					theseSkippedSNPs = append(theseSkippedSNPs, SNPNames[index])
					continue
				}
				correspondingSNP := SNPNames[index]
				hashAlleles[correspondingSNP] = append(hashAlleles[correspondingSNP], element)
			}
		}
		lineCounter += 1
	}
	return
}

//##############################################################################################################################

// This is the more complicated version of dat-parsing using file.Read() and a buffer
func parseDat(f string, hashAlleles map[string][]string, chromosome string, allSkippedSNPs map[string][]string) {
	var (
		file        *os.File
		properLine  string
		lineCounter int
		SNPNames    []string
		err         error
	)

	indexesOfSkippedSNPs := make(map[int]string) // Since go doesn't have the "set"-type, we have to cheat slightly

	if file, err = os.Open(f); err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	skippedSNPs := allSkippedSNPs[chromosome]

	// We can't use bufio.Scanner here, that leads to 'bufio.Scanner: token too long'
	buf := make([]byte, 1024)
	for {
		// read a chunk
		n, err := file.Read(buf)
		piece := string(buf[:n])
		// Is our current piece bridging a line? Break up the two lines and analyze the previous (now finished) line
		if strings.ContainsAny(piece, "\n") {
			twoPieces := strings.Split(piece, "\n")
			if len(twoPieces) != 2 {
				// The line is so short that we get several alleles per line. Shouldn't happen anymore?
				log.Fatal("FATAL ERROR: Tried to split up line into two, got ", len(twoPieces), " instead!\nFile is "+f+"\n")
			}
			// Add the rest of the previous line
			properLine += twoPieces[0]
			ll := strings.Split(strings.ToUpper(properLine), "\t")
			if lineCounter == 0 {
				// first line, check whether we have the correct file
				if ll[0] != "# FJFILE = GENOTYPE" {
					log.Fatal("FATAL ERROR: Parsing dat-file in dat-control-file, the file in question '" + f + "' doesn't seem to be a Flapjack file.\nHeader should be '# FJFILE = GENOTYPE'.")
				}
			} else if lineCounter == 1 {
				// second line, initialize the snps map
				SNPNames = ll[1:] // In flapjack-dat-headers, the first element is empty

				// Rename SNPNames form M1 to tapidor.c09|||M1
				for i, element := range SNPNames {
					SNPNames[i] = chromosome + SNPNameDelimiter + element
				}

				// We have to get the indexes for the SNPs in SNPNames so we can skip these in later lines
				for index, element := range SNPNames {
					// Is this element in our list of snps we should skip? If yes, don't include in positions array
					include := true
					for _, toCheck := range skippedSNPs {
						if toCheck == element {
							indexesOfSkippedSNPs[index] = "y"
							include = false
							break
						}
					}
					if include {
						_, ok := hashAlleles[element]
						if ok {
							log.Fatal("ERROR: The SNP ", element, " exists more than once in your dat-file! Aborting.")
						}
						hashAlleles[element] = make([]string, 0)
					}
				}
			} else {
				// All other lines - populate the hash with alleles
				// Don't need to keep the name of the individuals
				ll = ll[1:]
				// In some files, the lines are one longer than the header, with an empty element - try to delete this one
				if len(ll) != len(SNPNames) {
					ll = ll[:len(ll)-1]
				}

				// Still weird line - stop.
				if len(ll) != len(SNPNames) {
					log.Fatal("ERROR: The length of the line of SNP-names (line 2) is different than a line of alleles, SNP names are ", len(SNPNames), " length of alleles are in ", len(ll), " line number", lineCounter, "\nFile is ", f)
				}

				theseSkippedSNPs := make([]string, 0)
				for index, element := range ll {
					if isIntInMap(index, indexesOfSkippedSNPs) {
						theseSkippedSNPs = append(theseSkippedSNPs, SNPNames[index])
						continue
					}
					correspondingSNP := SNPNames[index]
					hashAlleles[correspondingSNP] = append(hashAlleles[correspondingSNP], element)
				}
			}
			// Start the next line
			properLine = twoPieces[1]
			lineCounter += 1
		} else {
			properLine += piece
		}
		if err != nil && err != io.EOF {
			log.Fatal(err)
		}
		// Are we done?
		if n == 0 {
			break
		}
	}
}

//##############################################################################################################################

func parseMap(f string) (skippedSNPs []string, nameSlice []string, posSlice []int) {
	var (
		file       *os.File
		chromosome string
		prevPos    int = -10000
		err        error
	)

	if file, err = os.Open(f); err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	scanner.Scan()
	if scanner.Text() != "# fjFile = MAP" {
		log.Fatal("FATAL ERROR: The Flapjack-map file '" + f + "' doesn't seem to be a map-file, header missing or wrong.\nShould be '# fjFile = MAP'.")
	}
	for scanner.Scan() {
		line := strings.Split(scanner.Text(), "\t")
		// Initialize the name of the chromosome
		if chromosome == "" {
			chromosome = line[1]
			// Test whether the output file, which is based on the name of the chromosome, exists
			f := chromosome + OutputSuffix
			_, err = os.Stat(f)
			if !os.IsNotExist(err) {
				log.Println("WARNING: '" + f + "' exists. Overwriting.")
			}
		}
		position, err := strconv.Atoi(strings.Replace(line[2], ",", "", -1))
		if err != nil {
			log.Fatal(err)
		}

		name := line[0]
		// Add the name of the chromosome for uniqueness
		name = chromosome + SNPNameDelimiter + name // no-one in the world uses ||| as a delimiter. I am looking forward to see this being broken.

		// We don't need to store SNPs <= X nt close to previous SNP, behave 99% identical anyway
		// SNPMinimumRange is a constant, see header of this file
		if position-SNPMinimumRange <= prevPos {
			skippedSNPs = append(skippedSNPs, name)
			continue
		}
		nameSlice = append(nameSlice, name)
		posSlice = append(posSlice, position)
		prevPos = position
	}
	return
}

//##############################################################################################################################

func placeContigs(unplacedMapDict map[string]*ContigMap, chromMapDict map[string]*ContigMap, unplacedHashAlleles map[string][]string, placedHashAlleles map[string][]string, comparisonMethod *string, cutoff *float64) {
	// We need a dict that shows which contig has been linked to what other contig

	var (
		scores_fh, unplaceableFh            *os.File
		err                                 error
		bestScore                           float64
		firstBestPartner, secondBestPartner *Contig
	)

	// Write table with stats for each placed contig
	log.Println("Writing scores to '" + ScoresFile + "'.")
	scores_fh, err = os.Create(ScoresFile)
	if err != nil {
		log.Fatal(err)
	}
	scores_fh.WriteString("Name of contig\tName of best partner\tScore\tSNPs\tNormalised recombinations\tStatus\n")
	defer scores_fh.Close()

	// Contigs that have no SNPs are unplaceable
	f := UnplaceableFile
	log.Println("Writing unplaceable contigs (no SNPs) to '" + f + "'.")

	if unplaceableFh, err = os.Create(f); err != nil {
		log.Fatal(err)
	}
	defer unplaceableFh.Close()

	unplaceableFh.WriteString("Name of contig\tReason\tNormalised recombinations\tNumber SNPs\n")

	for _, unplacedContigList := range unplacedMapDict {
		currentContig := unplacedContigList.firstContig
		// FormatFloat has format f (no exponent), 3 decimal points, float64

		for currentContig != nil {
			currentRecombinations := strconv.FormatFloat(currentContig.normalisedRecombinations, 'f', 3, 64)
			snps := currentContig.metaSNP
			numberOfSNPsOnContig := strconv.Itoa(len(currentContig.listOfSNPs))
			if len(snps) == 0 {
				// No SNPs on this contig
				unplaceableFh.WriteString(currentContig.name + "\tNo_SNPs\t" + currentRecombinations + "\t0\n")
				currentContig = currentContig.nextContig
				continue
			} else {
				// Find first and second best partner for the metaSNP
				// Iterate over all placed contigs
				firstBestPartner, bestScore, secondBestPartner = findBestPartner(snps, chromMapDict, comparisonMethod)

				// Write out best score
				scores_line := currentContig.name + "\t" + firstBestPartner.name + "\t"
				scores_line += strconv.FormatFloat(bestScore, 'f', 3, 64) + "\t" + numberOfSNPsOnContig + "\t" + currentRecombinations

				// With Hamming distance, we have to be below a cutoff, with all others, above
				if (*comparisonMethod == "hamming" && bestScore > *cutoff) || (*comparisonMethod != "hamming" && bestScore < *cutoff) {
					log.Println("REMOVING ", bestScore, *cutoff, currentContig.name)
					unplaceableFh.WriteString(currentContig.name + "\tScore_too_low\t" + currentRecombinations + "\t" + numberOfSNPsOnContig + "\n")
					currentContig = currentContig.nextContig
					scores_fh.WriteString(scores_line + "\tScore_too_low\n")
					continue
				}

				// Do we have a second best partner?
				if secondBestPartner != nil {
					// Two different best partners, into the trash it goes
					if secondBestPartner.chromosome != firstBestPartner.chromosome {
						unplaceableFh.WriteString(currentContig.name + "\tDifferent_partners\t" + currentRecombinations + "\t" + numberOfSNPsOnContig + "\n")
						currentContig = currentContig.nextContig
						scores_fh.WriteString(scores_line + "\tDifferent_partners\n")
						continue
					}
				}
				scores_fh.WriteString(scores_line + "\tPlaced\n")
				// We can place it!
				// Initialize the partner's map if needed
				if firstBestPartner.dictOfPartners == nil {
					firstBestPartner.dictOfPartners = make(map[float64][]*Contig)
				}

				// There is a chance, however rare, that two contigs have an identical score to the placed one. So we have to do slices
				partners, _ := firstBestPartner.dictOfPartners[bestScore]
				// Note: It doesn't matter whether the dict returns true or false, we always get a list (when empty, the dict also returns false)
				partners = append(partners, currentContig)
				firstBestPartner.dictOfPartners[bestScore] = partners

				// If the contig has more than 10 SNPs, there is a possibility that we can find an orientation using the last 1/3 of the contig
				if len(currentContig.lastMeta) > 0 {
					lastMeta := currentContig.lastMeta
					firstBestPartnerForLast, _, _ := findBestPartner(lastMeta, chromMapDict, comparisonMethod)
					if firstBestPartnerForLast != firstBestPartner {
						// Orientate
						// Rule: Does last_5_best_contig come before partnerFirstHalf, or does it come after?
						// If it comes before, the orientation is reverse complement
						firstStart := firstBestPartner.startPos
						last5Start := firstBestPartnerForLast.startPos
						if firstStart > last5Start {
							currentContig.orientation = "-"
						} else {
							currentContig.orientation = "+"
						}
					} else {
						// No orientation
						currentContig.orientation = "?"
					}
				} else {
					// No orientation
					currentContig.orientation = "?"
				}
			}
			// Go to next score
			currentContig = currentContig.nextContig
		}
	}
	return
}

//##############################################################################################################################

// Helper function - makes a Gff3 line from a Contig
func getGff3Line(contig *Contig) (lineSlice []string) {
	lineSlice = make([]string, 9)
	// Remember, gff3 is 1-inclusive
	start := strconv.Itoa(contig.startPos + 1)
	end := strconv.Itoa(contig.endPos + 1)
	lineSlice[0] = contig.chromosome
	lineSlice[1] = "fasta"
	lineSlice[2] = "contig"
	lineSlice[3] = start
	lineSlice[4] = end
	lineSlice[5] = "."
	lineSlice[6] = contig.orientation
	lineSlice[7] = "."
	lineSlice[8] = "ID=" + contig.name + ";Name=" + contig.name
	return lineSlice
}

//##############################################################################################################################

func writeGff3(chromMapDict map[string]*ContigMap, comparisonMethod *string) {

	for chrom, value := range chromMapDict {
		// Open a file-handle for this chromosome
		var (
			chromFh *os.File
			err     error
			f       string = chrom + OutputSuffix
		)
		log.Println("Writing gff3 for", f)

		// Try to create the file
		if chromFh, err = os.Create(f); err != nil {
			log.Fatal(err)
		}
		currentContig := value.firstContig

		chromFh.WriteString("##gff-version3\n")
		// TODO: Find a way to replace the -1. Not easy and not important.
		chromFh.WriteString("##sequence-region  " + currentContig.chromosome + "   1   -1\n")

		for currentContig != nil {
			lineSlice := getGff3Line(currentContig)
			chromFh.WriteString(strings.Join(lineSlice, "\t") + "\n")

			if len(currentContig.dictOfPartners) > 0 {
				// We can place unplaced contigs here!
				// sort the keys of the dictOfPartners by score
				scores := make([]float64, len(currentContig.dictOfPartners))
				i := 0
				for key, _ := range currentContig.dictOfPartners {
					scores[i] = key
					i++
				}
				// With Hamming, the lowest score is the best, so normal sorting is alright - in all other cases, high is better, so we have to invert
				if *comparisonMethod == "hamming" {
					// sort by lowest to highest
					sort.Sort(sort.Float64Slice(scores))
				} else {
					sort.Sort(sort.Reverse(sort.Float64Slice(scores)))
				}

				// Each value for the keys is a slice of contigs - we can't sort these, so keep them the way they showed up
				// In 99.9999% of cases the slice has only one contig
				for _, score := range scores {
					contigs := currentContig.dictOfPartners[score]
					for _, c := range contigs {
						lineSlice := getGff3Line(c)
						chromFh.WriteString(strings.Join(lineSlice, "\t") + "\n")
					}
				}
			}
			currentContig = currentContig.nextContig
		}
	}
}

//##############################################################################################################################

func main() {
	var name string
	// Get user arguments
	contigsFile := flag.String("contig", "", "The path to the file detailing contigs.")
	chromosomesFile := flag.String("chrom", "", "The path to the file detailing assembled chromosomes.")
	version := flag.Bool("v", false, "Prints the current version.")
	comparisonMethod := flag.String("comparison", "hamming", "Optional: Changes the SNP comparison mode. Options are 'hamming' (The standard. Penalized Hamming distance), 'percentage' (just normal comparison using percentage), 'r_squared' (r² from linkage disequilibrium).")
	cutoff := flag.Float64("cutoff", 0.0, "Optional: Cutoff for placement. Contigs with a Hamming score below this will not be placed.")

	flag.Parse()

	if *version {
		fmt.Println("Version:", CurrentVersion)
		os.Exit(0)
	}

	comparisonValid := isStrInMap(*comparisonMethod, ComparisonMethods)
	if !comparisonValid {
		log.Fatal("FATAL ERROR: Given comparison method is not a valid method (Valid methods are: 'hamming', 'percentage', 'r_squared'.")
	}

	if *cutoff == 0.0 && *comparisonMethod == "hamming" {
		*cutoff = math.MaxFloat64
	} // in all other cases it's 0

	if *contigsFile == "" || *chromosomesFile == "" {
		flag.Usage()
		os.Exit(1)
	}

	log.Println("Started.")
	log.Println("Contigs to be placed are in: '" + *contigsFile + "'.")
	log.Println("Already assembled chromosomes are in: '" + *chromosomesFile + "'.")

	// Parse map and dat file
	unplacedDats, unplacedMaps, unplacedGffs := parseControlFiles(*contigsFile)
	chromDats, chromMaps, chromGffs := parseControlFiles(*chromosomesFile)

	// Check gffs and maps whether there are any problems in naming of chromosomes
	checkForProblems(unplacedMaps, unplacedGffs)
	checkForProblems(chromMaps, chromGffs)
	checkDatFiles(unplacedDats)
	checkDatFiles(chromDats)

	// Parse the map-file for unplaced contigs
	// two slices, instead of one hash, since slices are ordered - less hassle to iterate over all of them
	// when we find out which SNPs are on which contig
	unplacedHashNames := make(map[string][]string)  // [tapidor.c09|||M1, tapidor.c09|||M2...]
	unplacedHashPositions := make(map[string][]int) // [10, 25...]
	allSkippedUnplacedSNPs := make(map[string][]string)
	listOfNames := make([]string, 0)
	log.Println("Now parsing unplaced maps.")
	for _, element := range unplacedMaps {
		skippedUnplacedSNPs, unplacedNames, unplacedPos := parseMap(element)
		name = strings.Split(unplacedNames[0], SNPNameDelimiter)[0] // tapidor.c09
		listOfNames = append(listOfNames, name)
		unplacedHashNames[name] = unplacedNames
		unplacedHashPositions[name] = unplacedPos
		allSkippedUnplacedSNPs[name] = skippedUnplacedSNPs
	}

	log.Println("Now parsing unplaced dats.")
	// Parse the dat-file for unplaced contigs
	unplacedHashAlleles := make(map[string][]string)
	for index, element := range unplacedDats {
		name = listOfNames[index]
		// Test whether the file is too long for go's "nice" parser -
		// If we get an error from bufio's Scanner, we have to use parseDat which uses a more complicated way to parse files
		// This is needed because bufio's Scanner has an upper-limit for lines, which is easily broken with .dats
		if linesAreTooLong(element) {
			parseDat(element, unplacedHashAlleles, name, allSkippedUnplacedSNPs)
		} else {
			parseDatBufio(element, unplacedHashAlleles, name, allSkippedUnplacedSNPs)
		}
	}
	// Do the same for placed contigs
	placedHashNames := make(map[string][]string)
	placedHashPositions := make(map[string][]int)
	allSkippedPlacedSNPs := make(map[string][]string)
	listOfNames = make([]string, 0)
	log.Println("Now parsing placed maps.")
	for _, element := range chromMaps {
		chromNames := make([]string, 0)
		skippedPlacedSNPs, chromNames, chromPos := parseMap(element)
		name = strings.Split(chromNames[0], SNPNameDelimiter)[0] // tapidor.C09|||M1 -> tapidor.C09
		listOfNames = append(listOfNames, name)
		placedHashNames[name] = chromNames
		placedHashPositions[name] = chromPos
		allSkippedPlacedSNPs[name] = skippedPlacedSNPs
	}
	// Parse the dat-file
	log.Println("Now parsing placed dats.")
	placedHashAlleles := make(map[string][]string) // Key: tapidor.c09|||M1, value:  [T,T,T,N,N,"",....]
	for index, element := range chromDats {
		name = listOfNames[index]
		if linesAreTooLong(element) {
			parseDat(element, placedHashAlleles, name, allSkippedPlacedSNPs)
		} else {
			parseDatBufio(element, placedHashAlleles, name, allSkippedPlacedSNPs)
		}
	}

	// Now parse the gff3-files, create a linked list of contigs for each chromosome
	log.Println("Now parsing gff3-files.")
	chromMapDict := make(map[string]*ContigMap)
	for _, element := range chromGffs {
		parseGffControlFile(element, placedHashPositions, placedHashNames, chromMapDict, placedHashAlleles)
	}
	unplacedMapDict := make(map[string]*ContigMap)
	for _, element := range unplacedGffs {
		parseGffControlFile(element, unplacedHashPositions, unplacedHashNames, unplacedMapDict, unplacedHashAlleles)
	}

	// Now iterate over the unplaced contigs, then iterate over all placed contigs and find the best partner
	log.Println("Now comparing all contigs.")
	placeContigs(unplacedMapDict, chromMapDict, unplacedHashAlleles, placedHashAlleles, comparisonMethod, cutoff)
	// Now iterate over all placed contigs, write them to gff3, and insert the unplaced contigs
	log.Println("Now writing GFF3-files.")
	writeGff3(chromMapDict, comparisonMethod)

	// Done!
	log.Println("Exit, pursued by a bear.")
	os.Exit(0)
}
