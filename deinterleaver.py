#! /bin/bash/env python
# Code to deinterleave FASTA files.
# Usage: python deinterleaver.py input_file output_file
import sys
linenum = 0
infile = open(sys.argv[1], 'r')
with open(sys.argv[2], 'w') as outfile:
    for line in infile:
        line = line.strip()
        if len(line) > 0:
            if line[0] == '>':
                if linenum == 0:
                    outfile.write(line + 'n')
                    linenum += 1
                else:
                    outfile.write('n' + line + 'n')
            else:
                outfile.write(line)
