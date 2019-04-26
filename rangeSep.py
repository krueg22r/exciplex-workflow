#! /usr/bin/env python3

#A program to determine the optimal range separation parameter for a given
#geometry. Operates on every .xyz file in the directory where it's run. 

import os
import glob
from scriptFunc import writeScriptLong

# list of range separation parameters to consider 
muList = [0.15, 0.20, 0.25, 0.30, 0.35]

structList = glob.glob('*xyz')
for struct in structList:		      
	sysName = struct.replace('.xyz', '') 
	#print shell script to submit jobs at one time 
	submitFile = open(sysName+'_submit.sh', 'w')
	print('#! /bin/bash', file=submitFile)
	os.system('mkdir '+sysName)
	#write orca input file
	infile = open('mon', 'w')
	print('! LC-BLYP def2-TZVP TightSCF Grid5 NoFinalGrid', file=infile)
	print('', file=infile)
	print('%MaxCore 1800', file=infile)
	print('', file=infile) 
	print('%pal nprocs 12', file=infile)
	print('end', file=infile) 
	print('', file=infile) 
	print('%method', file=infile)
	print('	RangeSepMu num', file=infile)
	print('end', file=infile)
	print('', file=infile)
	print('* xyz charge spin ', file=infile) #charge and spin later replaced with real vals
	#print initial coordinates from xyz file
	coordFile = open(struct, 'r') 
	coordFile.readline()
	coordFile.readline()
	for line in coordFile: 
		print(line, file=infile, end='')
	print('*', file=infile) 
	infile.close()
	# write queue script. give base name, use default 12 procs 
	writeScriptLong('mon') 
	for mu in muList: 
		os.system('mkdir '+sysName+'/'+str(mu))
		chargeList = ['0', '1']
		for charge in chargeList:
			dirName = sysName+'/'+str(mu)+'/'+charge
			os.system('mkdir '+dirName)
			# move files into directory 
			os.system('cp mon.sh mon '+dirName)
			fileName = dirName+'/mon'
                        # edit files to make spin and charge accurate 
			os.system('sed \'s/num/'+str(mu)+'/\' --in-place '+fileName)
			spin = 5
			if int(charge) == 0:
				spin = 1
			else:
				spin = 2
			os.system('sed \'s/spin/'+str(spin)+'/\' --in-place '+fileName)
			os.system('sed \'s/charge/'+str(charge)+'/\' --in-place '+fileName)
			# write details to job submission file
			print('cd '+dirName, file = submitFile)
			print('qsub mon.sh', file = submitFile)
			print('sleep 2s', file=submitFile) 
			print('cd ../../..', file = submitFile) 
        # remove files from top directory once they are in place 
	os.system('rm mon mon.sh')
	submitFile.close()
