#! /usr/bin/env python
# With two files, one containing concatenated sequences and the other gene ranges,
# separate all into their own sequences in FASTA format.

# Presumably copying and pasting from a NEXUS file ...
# Example of sequence_file:
# Taxon1 MPPLS---SLDG ...
# Taxon2 MPPLSVPLSLDE ...
# Taxon3 MPSL--PLSLDG ...

# Example of generange_file:
# cytb:1-2i9,
# COX1:220-608,
# NAD5:609-1051,
# ...

# Usage: python deconcatenator.py sequence_file generange_file output_name
import sys

genes_and_ranges = {}
with open(sys.argv[2], 'r') as genefile:
    for line in genefile:
        line = line.strip().strip(',').split(':')
        line[1] = line[1].split('-')
        genes_and_ranges[line[0]] = [line[1][0], line[1][1]]

sequences = {}
with open(sys.argv[1], 'r') as myxofile:
    for line in myxofile:
        line = line.strip().split('\t')
        for gene in genes_and_ranges:
             sequences['>' + line[0] + '_' + gene] = line[1][(int(genes_and_ranges[gene][0])-1):int(genes_and_ranges[gene][1])]

with open(sys.argv[3], 'w') as outfile:
    for i in sequences:
        outfile.write(str(i) + '\n' + str(sequences[i].replace('-','') + '\n'))
#       If you don't want to unalign sequences:
#       outfile.write(str(i) + '\n' + str(sequences[i] + '\n')) 
