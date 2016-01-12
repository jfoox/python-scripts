#! /usr/bin/env python
import urllib
import os
import sys
# This a script to calculate percent non-coding regions
# in mitochondrial genomes: a comically difficult task for 
# something seemingly so simple.
# USAGE:
# python igr_calculator.py accessionsfile

# First step is to have a file of accession numbers,
# one on each line, and turn that into a list.
accessions = []
with open(sys.argv[1], 'r') as acc_file:
	for line in acc_file:
		line = line.strip()
		accessions.append(line)
	acc_file.close()

# Now, for each accession, download its respective .gb file from NCBI.
for accession in accessions:
	url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=gb&id=' + accession
	print 'Downloading ' + url
	urllib.urlretrieve(url, 'sequences/' + accession + '.gb')

### The next step is to convert each .gb file into a .bed file.
### Unfortunately, this appears to be only doable manually through a program like Geneious.
### Once you have respective .bed files, the fun begins.

# Take each .bed file, find only certain lines, write those to temporary files.
for accession in accessions:
	lines = open('sequences/' + accession + '.bed').read().splitlines()
	for line in lines:
		with open('sequences/' + accession + '.tmp', 'a') as tempfile:
			if 'CDS' in line or 'intron' in line or 'RNA' in line or 'source' in line:
				tempfile.write(line + '\n')
		tempfile.close()

# Now we have temporary files with just lines of interest.
# Next, we save the taxon name as a key,
# then do all the calculations, and save the sum as the key's value.
for accession in accessions:
	with open('IGR_calculations.txt', 'a') as endfile:
		lines = open('sequences/' + accession + '.tmp').read().splitlines()

# Start with a list of numbers the size of the genome (basically enumerate nucleotides),
# with the null assumption that they are all non-coding,
# then selectively remove coding regions based on .bed files.

# First, figure out how big the genome is.
		for line in lines:
			line = line.split('\t')
			if 'source' in line[3]:
				totalsize = line[2]
		genome = range(0, int(totalsize) + 1)
		
# Now we have a big list of nucleotides that we will iteratively remove.
# First, we'll deal with introns, because they suck the hardest.

# We want to find where the intron starts and ends. Then find ranges in between.
# Then: (lowest_inbetween:intron_start), (intron_end:highest_inbetween)
# will give us those lovely flanking intron regions to pluck out.
		for line in lines:
			line = line.strip().split('\t')
			tweeners = []
			if 'intron' in line[3]:
				intron_first = int(line[1])
				intron_last = int(line[2])
				for line in lines:
					line = line.strip().split('\t')
					if intron_last > int(line[1]) > intron_first:
						tweeners.append(int(line[1]))
						tweeners.append(int(line[2]))
				tweeners.sort()
				print tweeners
				if len(tweeners) > 0:					
					intron_left = range(intron_first, tweeners[0])
					intron_right = range(tweeners[-1], intron_last)
					for i in intron_left:
						genome[i] = '!'
					for i in intron_right:
						genome[i] = '!'

# But sometimes, the "intron" is not labeled, so we have to check ND5 itself...
		ND5ers = []
		ND5tweeners = []
		for line in lines:
			if 'intron' not in lines:
				line = line.strip().split('\t')
				if 'ND5' in line[3]:
					ND5ers.append(int(line[1]))
					ND5ers.append(int(line[2]))
					ND5ers.sort()
		print 'ND5 ranges: ' + str(ND5ers)
		for line in lines:
			line = line.strip().split('\t')
			if 'ND5' not in line[3] and len(ND5ers) > 0 and ND5ers[-1] > int(line[1]) > ND5ers[0]:
				ND5tweeners.append(int(line[1]))
				ND5tweeners.append(int(line[2]))
		ND5tweeners.sort()
		print 'Stuff in between ND5: ' + str(ND5tweeners)
		if len(ND5tweeners) > 0: 
			intron_left = range(ND5ers[0], ND5tweeners[0])
			intron_right = range(ND5tweeners[-1], ND5ers[-1])
			for i in intron_left:
				genome[i] = '!'
			for i in intron_right:
				genome[i] = '!'
				
# Now we will take care of all other coding regions.
		for line in lines:						
			line = line.split('\t')
			if 'source' in line[3]:
				endfile.write(line[3] + '\t')
				print line[3]
			elif 'CDS' in line[3]:
    			# Need to add 1 to first value because for some reason gb->bed subtracts a number
				CDS = range(int(line[1]) + 1, int(line[2]) + 1)
				for i in CDS:
					genome[i] = '!'
			elif 'RNA' in line[3]:
				RNA = range(int(line[1]) + 1, int(line[2]) + 1)
				for i in RNA:
					genome[i] = '!'
					
# Now that all coding regions are accounted for, time to calculate IGR.
		if genome[0] != '!':
			genome.remove(0)
		for i in range(0, genome.count('!')):
			genome.remove('!')
		IGR = '%.2f' % (float(len(genome)/float(totalsize))*100)
		print 'The amount of noncoding for accession # %s is %s%%' % (accession, IGR)
		endfile.write(IGR + '\n')
