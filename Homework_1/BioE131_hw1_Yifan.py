#! /usr/bin/env python3
####Plz run this script in python 3
##################
####Reorganize the FASTA files
####Yifan Ethan Chen
####BioE 131 HW1
####09/20/17 Cal


import sys
from sys import exit
from sys import argv

def conseq(sequence):
    
    sequence = sequence.replace('A', '1')
    sequence = sequence.replace('a', '2')
    sequence = sequence.replace('T', 'A')
    sequence = sequence.replace('t', 'a')
    sequence = sequence.replace('U', 'A')
    sequence = sequence.replace('u', 'a')
    sequence = sequence.replace('1', 'T')
    sequence = sequence.replace('2', 't')
    
    sequence = sequence.replace('C', '1')
    sequence = sequence.replace('c', '2')
    sequence = sequence.replace('G', 'C')
    sequence = sequence.replace('g', 'c')
    sequence = sequence.replace('1', 'G')
    sequence = sequence.replace('2', 'g')
    
    sequence = sequence.replace('R', '1')
    sequence = sequence.replace('r', '2')
    sequence = sequence.replace('Y', 'R')
    sequence = sequence.replace('y', 'r')
    sequence = sequence.replace('1', 'Y')
    sequence = sequence.replace('2', 'y')
    
    sequence = sequence.replace('K', '1')
    sequence = sequence.replace('k', '2')
    sequence = sequence.replace('M', 'K')
    sequence = sequence.replace('m', 'k')
    sequence = sequence.replace('1', 'M')
    sequence = sequence.replace('2', 'm')
    
    sequence = sequence.replace('B', '1')
    sequence = sequence.replace('b', '2')
    sequence = sequence.replace('V', 'B')
    sequence = sequence.replace('v', 'b')
    sequence = sequence.replace('1', 'V')
    sequence = sequence.replace('2', 'v')
    
    sequence = sequence.replace('D', '1')
    sequence = sequence.replace('d', '2')
    sequence = sequence.replace('H', 'D')
    sequence = sequence.replace('h', 'd')
    sequence = sequence.replace('1', 'H')
    sequence = sequence.replace('2', 'd')
    
    sequence = sequence[::-1] 
    return sequence

def writeout(proteinname,sequence):
    fout.write(proteinname+'\n')
    seqlist = [sequence[i:i+80]+'\n' for i in range(0,len(sequence),80)]
    for i in seqlist:
        fout.write(i)

###Error Handle
try: 
    filename = sys.argv[1]
except IndexError:
    print ('The user did not specify a filename argument on the command line.')
    sys.exit()
try:
    fin = open(filename)
except IOError:
    print ('The file named by the user does not exist, or is not readable.')
    sys.exit()
    
contents = fin.readlines()
fout = open(filename+'_output.fa', 'w')

n = 0
while n < len(contents):
    if contents[n] == '\n':
        n += 1
        pass
    elif contents[n][0] == '>':
        proteinname = contents[n].split()[0]
        sequence = ''
        n += 1
        while contents[n][0] == '\n':
            n += 1
        if contents[n][0] == '>':
            print ('No sequence for ' + proteinname[1:] + '.')    
        while n < len(contents) and contents[n][0] !='>':
            if contents[n][-1:] != '\n':
                sequence += contents[n][0:]
            else:
                sequence += contents[n][0:-1]
            n += 1
        for i in sequence:
            if i in ['E', 'F', 'I', 'L', 'P', 'Q', 'e', 'f', 'i', 'l', 'p', 'q']:
                print ('The file named by the user contains sequences that are not DNA.')
                sys.exit()
        sequence = conseq(sequence)
        writeout(proteinname,sequence)
    elif contents[n][0] != '>':
        print ('The file named by the user is not in FASTA format.')
        sys.exit()

fout.close()
sys.exit()