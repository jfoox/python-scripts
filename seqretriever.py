#! /usr/bin/env python
import sys
import os
# A script for extracting a certain subset of sequences from within a FASTA file.
# Usage: python seqretriever.py fasta_file sequence_list_file
 
# First, convert FASTA file into file with one line per sequence.
# Make sure the name of your FASTA filename doesn't contain any dots besides the one before the extension!
firstline = True
onelinefile = 'allononeline_' + sys.argv[1].partition('.')[0] + '.tmp'
with open(onelinefile, 'w') as outfile:
    for line in open(sys.argv[1], 'r'):
        line = line.strip()
        if len(line) > 0:
            if line[0] == '>':
                if firstline == True:
                    outfile.write(line + '\t')
                    firstline = False
                else:
                    outfile.write('\n' + line + '\t')
            else:
                outfile.write(line)
 
# Populate a dictionary using the 'allononeline' file
all_seqs = {}
with open(onelinefile, 'r') as allseqsfile:
    for line in allseqsfile:
        line = line.strip().split('\t')
        all_seqs[line[0][1:]] = line[1]
 
# Generate a set of the sequences you wish to retrieve
desired_seqs = set()
with open(sys.argv[2], 'r') as listfile:
    for line in listfile:
        desired_seqs.add(line.strip())
 
# Find the overlap between the total sequences and the desired ones
all_seqs_names = set(all_seqs.keys())
toextract = all_seqs_names.intersection(desired_seqs)
 
# Use 'toextract' set to generate desired file
with open('justdesired_' + sys.argv[1], 'w') as extractfile:
    for name in toextract:
        extractfile.write('>' + name + '\n' + all_seqs[name] + '\n')

os.remove(onelinefile)
