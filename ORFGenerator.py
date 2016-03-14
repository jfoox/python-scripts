#! /usr/bin/env python
# ORF Generator
# python ORFG.py -in input.fa -out output.fa -deint y (default=y) -length 
import os
import sys
import argparse

# command line input
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='ORF Generator. Fun!')
parser.add_argument('-i', '--input', metavar='', default='', help='Input FASTA file; specify path if not in same folder')
parser.add_argument('-o', '--output', metavar='', default='', help='Output file name')
parser.add_argument('-l', '--length', metavar='', type=int, default=20, help='Minimum length of translated amino acid sequence to retain (default: 20)')
args = parser.parse_args()


codondict = {'ATT':'I',	'ATC':'I',	'ATA':'I',	'CTT':'L',	'CTC':'L',	
'CTA':'L',	'CTG':'L',	'TTA':'L',	'TTG':'L',	'GTT':'V',	'GTC':'V',	
'GTA':'V',	'GTG':'V',	'TTT':'F',	'TTC':'F',	'ATG':'M',	'TGT':'C',	
'TGC':'C',	'GCT':'A',	'GCC':'A',	'GCA':'A',	'GCG':'A',	'GGT':'G',	
'GGC':'G',	'GGA':'G',	'GGG':'G',	'CCT':'P',	'CCC':'P',	'CCA':'P',	
'CCG':'P',	'ACT':'T',	'ACC':'T',	'ACA':'T',	'ACG':'T',	'TCT':'S',	
'TCC':'S',	'TCA':'S',	'TCG':'S',	'AGT':'S',	'AGC':'S',	'TAT':'Y',	
'TAC':'Y',	'TGG':'W',	'CAA':'Q',	'CAG':'Q',	'AAT':'N',	'AAC':'N',	
'CAT':'H',	'CAC':'H',	'GAA':'E',	'GAG':'E',	'GAT':'D',	'GAC':'D',	
'AAA':'K',	'AAG':'K',	'CGT':'R',	'CGC':'R',	'CGA':'R',	'CGG':'R',	
'AGA':'R',	'AGG':'R',	'TAA':'*',	'TAG':'*',	'TGA':'*'}


def deinterleaver(filename):
    '''
    Deinterleave FASTA file if not already done so. Necessary for parsing through file.
    '''
    linenum = 0
    with open('ORFGtemp_' + filename, 'w') as outfile:
        with open(filename, 'r') as infile:
            for line in infile:
                line = line.strip()
                if len(line) > 0:
                    if line[0] == '>':
                        if linenum == 0:
                            outfile.write(line + '\n')
                            linenum += 1
                        else:
                            outfile.write('\n' + line + '\n')
                    else:
                        outfile.write(line)

def translator(current_seq, current_seq_name, frame):
    '''
    For the given reading frame, translate the nucleotide sequence into amino acids.
    '''
    translated = []
    if frame < 0:
        current_seq = current_seq[::-1]
        frame = abs(frame)
    for i in range(frame - 1, len(current_seq), 3):
        codon = current_seq[i:i+3]
        if len(codon) == 3:
            if codon in codondict:
                translated.append(codondict[codon])
            else:
                sys.exit('Error: invalid codon \"%s\" in sequence \"%s\"' % (codon, current_seq_name))
    return ''.join(translated)
    

# Check to ensure that the FASTA file is deinterleaved. If not, call deinterleaver function.
file_to_use = args.input
with open(args.input) as checkfile:
    top = [next(checkfile) for x in xrange(3)]
if top[0][0] == '>' and top[2][0] == '>':
    pass
else:
    deinterleaver(args.input)
    file_to_use = args.input + 'ORFG.tmp'


frames = [1, 2, 3, -1, -2, -3]
with open(args.output, 'w') as outfile:
    for line in open(file_to_use, 'r'):
        each_frame_trans = {}
        line = line.strip()
        if line[0] == '>':
            current_seq_name = line[1:]
        if line[0] != '>':
            for frame in frames:
                translated_at_that_frame = translator(line, current_seq_name, frame)
                each_frame_trans[frame] = translated_at_that_frame
    
            # Now look at translated sequence and determine which part is ORF.
            longest_chunks = {}
            chunk_comp = []
            for translation in each_frame_trans:
                longest_chunks[sorted(each_frame_trans[translation].split('*'), key=len)[-1]] = translation
            for chunk in longest_chunks:
                chunk_comp.append(chunk)
            best_chunk = sorted(chunk_comp, key=len)[-1]

            if len(best_chunk) >= args.length:
                outfile.write('>' + current_seq_name + '_Frame_' + str(longest_chunks[best_chunk]) + '\n' + best_chunk + '\n')


if file_to_use.endswith('ORFG.tmp'):
    os.remove(file_to_use)
print 'Completed! Sequences written to: %s' % args.output
