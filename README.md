# contigPlacer

It places contigs!

## What does it do?

It compares unplaced contigs to SNPs/genotypes on placed contigs and places unplaced contigs next to their best partner using a penalised Hamming distance. In a way, it's similar to genetic mapping algorithms - but genetic mapping algorithms make de novo maps, I just wanted to place some newly assembled contigs into an already existing assembly.

It takes genotyping data as input.

### What does it do, in detail?

For each unplaced contig, it gets all SNPs located on that SNPs and merges them into one "metaSNP" which for each individual carries the most common allele. When using genotypes, recombinations on a contig are rare, so in many cases the genotypes for one individual on one contig are identical anyway. It also merges all SNPs for placed contigs.

Then, it compares all unplaced contigs against all placed contigs with a penalised Hamming distance - if both alleles for the individual are missing, add 0.75, if one is missing, add 0.5, if they are different, add 1, if they are identical, add 0. That means that the contig placed 

Then it does a few extra tricks to increase accuracy:

1. Are the best and the second-best placed partner contig on different chromosomes? If so, discard the contig, it's not placeable.

2. Can we orient the contig? Take the first and last 5 SNPs if there are more or equal to 10 SNPs located on it, and see whether the last 5 SNPs are more similar to the SNPs of the preceding contig, and whether the first 5 SNPs are closer to the succeeding contig. If that is the case, reverse complement the contig. (Note: This happens very rarely).

## How long does it run?

With a large plant genome and about a million SNPs it takes a few hours using one CPU on my laptop.

## How accurate is it?

Probably not *that* accurate - some (unpublished) tests show that the chromosome placement is ~ 98% correct, but I cannot test the base-pair level placement without optical maps etc. I would assume that the base-pair level placement is not very accurate.

## What does it take as input files?

For each collection of unplaced contigs, and for each pseudo-molecule assembly, contigPlacer takes three input files.

Since I work mostly with [Flapjack](https://ics.hutton.ac.uk/flapjack/), I used the input map and dat file from Flapjack, and a gff3 file detailing where each contig is. A simple example:

The map file (tab-delimited), one SNP per line, SNP-ID, chromosome, position:

```
#fjFile = MAP
M1  A01 10
M2  A01 20
M1000 A01   9,292
etc.
```

The dat file (tab-delimited), one individual per line, first the header with all SNP-IDs, then for each line individual-ID, and all genotypes:

```
#fjFile = GENOTYPE
    M1  M2  M1000 etc.
Ind1    A    A    B etc.
Ind2    B    B    A etc.
```

The gff3 file detailing contigs, one contig per line, only the 4th, 5th and 9th element are important:

```
##gff-version    3
##sequence-region        A01    1       29136790
A01 fasta   contig  1   4521    .   .   .   ID=Contig_1;Name=Contig_1
A01 fasta   contig  4522    8999    .   .   .   ID=Contig_2;Name=Contig_2
etc.
```

The paths to these control files are detailed in two control files - one file for the unplaced contigs, and one file for the placed contigs in pseudo-chromosomes, always in the order dat, map and gff3 (tab-delimited)

Placed contigs file:

```
placed_A01.dat placed_A01.map placed_A01.gff3
placed_A02.dat placed_A02.map placed_A02.gff3
etc.
```

And the same for the unplaced contigs file. Relative or absolute paths are fine.

## Dependencies

Install Go from https://golang.org/

I've tested it with a few versions from 1.2 to 1.5.1. I've tried to keep it free from external dependencies as that just makes things complicated.

I've only tested it under various Linux distributions, no clue whether it works under Windows.

## Installation

Run `go build` in this directory to get the compiled contigPlacer. As it's a self-contained binary you can deploy the resulting file to your server.

## Usage

Either

    go run contigPlacer.Go

or (after building)

    ./contigPlacer

with these flags:

    -chrom chromFile -contig toPlaceFile

Both `chromFile` and `toPlaceFile` are files detailing paths to the placed and unplaced contigs - see above.

Optional flags:

    -cutoff some_number 

This does not place unplaced contigs when their best partner score is above a certain threshold - I'd use 1/3 * the size of the population you have. I STRONGLY recommend you to use this flag.

    -comparison hamming (default), percentage, r_squared

This changes the comparison function. By default, it uses the Hamming distance. A few tests show this option not to make a big difference - since Hamming is penalising missing alleles I'd assume you place less in unimputed datasets compared to simple `percentage`.

## Output

You get one gff3 file per collection of placed contigs detailing where the original contigs are located, and where the new contigs are placed. You can use `makeChromosomesFromContigplacer.py` to make fasta files out of these files.

You also get one file called `all_scores.txt` which details the best score (Hamming distance - 0 means no difference between the unplaced contig and the best partner-contig) for each unplaced contig. This file also stores for each unplaced contig how many recombinations per individual were counted on the contig.

You also get `list_of_unplaceable_contigs.txt`, which details for what reason which contig couldn't be placed. This can be "No\_SNPs", "Different\_partners" (first and second best partner of an unplaced contig are on different chromosomes), and Score\_too\_low when you've supplied a cutoff.
