# -*- coding: utf-8 -*-
''' 
This python code searches for two very conserved sequences within the alignmnet
and returns an integer of the position where the aa of interest is

Daniel Martinez, August - 2021
'''

from Bio import SeqIO
import sys

# input and output
inFile = sys.argv[1]

# read fasta file
records = list(SeqIO.parse(inFile, "fasta"))

### functions
def pattern_find(seqs):
	temp, pos1, pos2 = 0, [], []
	# find patterns in sequences
	for i in range(len(records)):
		if ('GVEA' in records[i].seq)  == True:
			temp = temp + 1
			pos1.append(records[i].seq.find('GVEA') + 4)
		elif ('GIEA' in records[i].seq) == True:
			temp = temp + 1
			pos2.append(records[i].seq.find('GIEA') + 4)
		else:
			pass
	if temp > 0:
		print('Patterns found in ', str(temp), 'sequences!\n')
	elif temp == 0:
		print('Pattern not found...\n')
		break
	else:
		print('Something went wrong.\n')
	print(('For GVEA pattern, the aa is at position: ' + str(list(set(pos1))[0] + 1)))
	print(('For GIEA pattern, the aa is at position: ' + str(list(set(pos2))[0] + 1)+ '\n'))
	pos1 = list(set(pos1))[0] + 1
	print(pos1)
	sys.exit(0)
	# return(pos1)

pattern_find(records)