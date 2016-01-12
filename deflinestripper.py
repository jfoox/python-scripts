#! /usr/bin/env python
# The output of Trinity generates enormous definiton line names
# that often confuse mapping algorithms like bwa.
# This strips long deflines but leaves the unique names.
# Usage: python deflinestripper.py filename

import sys
linenum = 0
with open('shortdeflines_' + sys.argv[1], 'w') as outfile:
	with open(sys.argv[1], 'r') as infile:
		for line in infile:
			line = line.strip()
			if len(line) > 0:
				if line[0] == '>':
					seqname = line.partition(' ')[0][1:]
					if linenum == 0:
						outfile.write('>' + sys.argv[2] + '_' + seqname + '\n')
						linenum += 1
					else:
						outfile.write('\n' + '>' + sys.argv[2] + '_' + seqname + '\n')
				else:
					outfile.write(line)
