#! /usr/bin/env python3

# after tuning is performed using rangeSep.py, this program 
# extracts the data from those calculations and sets up an 
# excited-state optimization calculation for each structure 
# in the directory 

import os
import glob
from muFuncs import getMu
from scriptFunc import writeScriptLong
import argparse

# open file for shell script to submit all jobs 
submitFile = open('batch_submit.sh', 'w') 
structList = glob.glob('*xyz')
parser = argparse.ArgumentParser()
# user has the option to provide mu value, rather than 
# extracting it from tuning runs in same dir 
parser.add_argument('-mu', '--mu', action = 'store_true')
args = parser.parse_args()
for struct in structList: 
	mu = 0
	dirName = struct.replace('.xyz', '')
	if args.mu: 
		mu = float(input('enter mu for struct '+struct+'\n')) 	
	else: 
		mu = getMu(dirName)
	optDir = dirName+'/S1opt'
	os.system('mkdir '+optDir)	
	# write input for Orca 
	infile = open('complex', 'w')
	print('! RKS LC-BLYP def2-TZVP TightSCF Grid5 NoFinalGrid Opt', file=infile)
	print('', file=infile)
	print('%MaxCore 1800', file=infile)
	print('', file=infile)
	print('%pal nprocs 12', file=infile)
	print('end', file=infile)
	print('', file=infile)
	print('%tddft', file=infile)
	print('NRoots 2', file=infile) 
	print('MaxDim 8', file=infile) 
	print('Iroot 1', file=infile) 
	print('end', file=infile)
	print('', file=infile)
	print('%method', file=infile)
	print(' RangeSepMu '+str(round(mu,4)), file=infile) 
	print('end', file=infile)
	print('', file=infile)
	if 'adical' in dirName: #if open-shell ground state
		print('*xyz 0 2', file=infile)
	else: 
		print('* xyz 0 1 ', file=infile)
	# print coordinates to input file 
	coordFile = open(struct, 'r')
	coordFile.readline()
	coordFile.readline()
	for line in coordFile:
		print(line, file=infile, end='')
	print('*', file=infile)
	infile.close()
	writeScriptLong('complex') #write queueing system infile 
	os.system('mv complex complex.sh '+optDir)
	print('cd '+optDir, file=submitFile) 
	# finish writing shell script 
	print('qsub complex.sh', file = submitFile) 
	print('sleep 2s', file=submitFile)
	print('cd ../..', file = submitFile)
submitFile.close()
	
	
