from Bio import AlignIO
import re
import sys
import argparse

informats = ["fas", "a2m", "a3m", "sto", "psi", "clu"]
outformats = ["fas", "a2m", "a3m", "sto", "psi", "clu", "ufas"]

formats_translate = {'fas':'fasta',
                     'a2m':'fasta',
                     'a3m':'fasta',
                     'sto':'stockholm',
                     'clu':'clustal',
                     'psi':'fasta',
                     'ufas':'fasta',
                     }

parser = argparse.ArgumentParser(description='Reformat a file')
parser.add_argument('informat', help='format of the input file', choices=informats, type=str) #required=True, 
parser.add_argument('outformat', help='format of the output file', choices=outformats, type=str) # required=True,
parser.add_argument('infile', help='input file', type=str)
parser.add_argument('outfile', help='output file', type=str)
args = parser.parse_args(sys.argv[1:])

informat = args.informat
outformat = args.outformat
infile = args.infile
outfile = args.outfile

# Assign informat, outformat, infile, and outfile
if not outformat:
	outformat_from_file = outfile.split('.')[-1].lower()
	if outformat_from_file in ['fa', 'fasta', 'afa', 'afas', 'afasta']:
		outformat = 'fas'
	elif outformat_from_file == 'aln':
		outformat = 'clu'
	else:
		print("Using FASTA output format: %s\n" %outformat_from_file)
		outformat="fas"

if not informat:
	informat_from_file = infile.split('.')[-1].lower()
	if informat_from_file == 'aln':
		informat = 'clu'
	elif informat_from_file in ['fa', 'fasta']:
		informat = 'fas'
	else:
		print("Using FASTA input format: %s\n" %informat_from_file)
		informat="fas"

################################################################################################
# Reformat a single file
################################################################################################
alignments = AlignIO.read(infile, formats_translate[informat])

################################################################################################
# Output part
################################################################################################
AlignIO.write(alignments, outfile, formats_translate[outformat])