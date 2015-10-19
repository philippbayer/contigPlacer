#!/usr/bin/env python
import sys
VERSION = 1.0
if "-v" in sys.argv:
    sys.stderr.write(VERSION+"\n")
    sys.exit(0)

try:
    fasta = sys.argv[1]
    files = sys.argv[2:]
except:
    sys.stderr.write("Creates new FASTA-files from the output of contigPlacer.\n")
    sys.stderr.write("-v prints version number\n")
    sys.stderr.write("First argument is a file containing ALL contigs, then come n arguments for all gff3-files that contigPlacer produced.\n")
    sys.stderr.write("Usage: python makeChromosomesFromContigplacer.py file_with_all_contigs.fasta first.gff3 second.gff3 third.gff3 ...\n")
    sys.stderr.write("or:    python makeChromosomesFromContigplacer.py file_with_all_contigs.fasta *.gff3\n")
    sys.exit(1)

from Bio import SeqIO

index = SeqIO.index(fasta, "fasta")
seen = set()

for el in files:
    fh = open(el)
    out = el.replace(".gff3", ".fasta")
    sys.stderr.write("Writing to %s\n"%out)
    out_f = open(out, "w")
    for line in fh:
        if "##" in line: continue
        ll = line.rstrip("\n\r").split("\t")
        orientation = ll[6]
        id = ll[-1].split(";")[0].replace("ID=","")
        seen.add(id)
        seq = index[id]
        if orientation == "-":
            # reverse complement
            seq.seq = seq.seq.reverse_complement()
        out_f.write(seq.format("fasta"))

unplaced_out = open("Unplaced_contigs.fasta", "w")
sys.stderr.write("Writing all unplaced contigs to Unplaced_contigs.fasta\n")
for element in index:
    if element not in seen:
        unplaced_out.write(index[element].format("fasta"))

