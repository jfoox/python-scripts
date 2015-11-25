#! /bin/bash/env python
# Code to deinterleave FASTA files.
# Usage: python deinterleaver.py file_to_deinterleave output_file_name
import sys
linenum = 0
infile = open(sys.argv[1], 'r')
with open(sys.argv[2], 'w') as deintfile:
    for line in infile:
        line = line.strip()
        if len(line) > 0:
            if line[0] == '>':
                if linenum == 0:
                    deintfile.write(line + 'n')
                    linenum += 1
                else:
                    deintfile.write('n' + line + 'n')
            else:
                deintfile.write(line)
