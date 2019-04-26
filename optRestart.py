#! /usr/bin/env python3

# Unfortunately, the excited state opt may not finish within the 
# wall clock limit. This program sets up a calculation using the 
# latest structure obtained from the optimization as init coordinates

import os
import glob
from scriptFunc import writeScriptLong

trajName = glob.glob('*trj') 
trajName = max(trajName, key=len) #in case it's been restarted before...
baseName = trajName.replace('.trj', '') 
#open necessary files
origInfile = open(baseName, 'r') 
newInfile = open(baseName+'Restart', 'w')
coordFile = open(trajName, 'r') 
nAtoms = int(coordFile.readline())
#collect last coordinates
os.system('tail -'+str(nAtoms)+' '+trajName+'  > tmp.xyz')
for line in origInfile: 
#copy the first section of the old infile that doesn't change
	if 'xyz' not in line: 
		print(line, file=newInfile, end='') 
	else: # xyz in line signals that we've reached the coord part 
		print(line, file=newInfile, end='')
		break
# make a temp file with final xyz coords from last run 
endCoords = open('tmp.xyz', 'r') 
#print these coords to the new infile 
for line in endCoords: 
	print(line, file=newInfile, end='')
print('*', file=newInfile) 
#write queue script
writeScriptLong(baseName+'Restart')	
