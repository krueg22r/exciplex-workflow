#! /panfs/panbox/home/rkrueger/python/bin/python3

# a script capable of matching the complexes w/data about their monomers 

import glob
import numpy as np
from itertools import combinations 

def printHeader(): 
	print('# column data | 0: complex name, 1: homodimer, 2: oxygenated, 3: complex heavy atom mass, 4: num overlapping heavy atoms')
	print('# 5: large mon planarity index; 6: small mon planarity index, 7: large mon mass, 8: small mon mass')
	print('# 9: large mon bandgap, 10: small mon bandgap, 11: large mon emission')
	print('# 12: small mon emission, 13: complex emission')
	print('# 14: diameter large mon, 15: diameter small mon')
	print('# 16: avg interatomic distance, 17: avg bond length, 18: La?, 19: Large mon Avg bond len, 20: Small mon avg bond len')
# get avg intermolecular distance 
getIntermolDist(coordsBig, coordsSmall): 
	avgIntermolDist = 0
	nPairs = 0
	for i in range(0, coordsBig.shape[0]):
		 for j in range(0, coordsSmall.shape[0]):
			dist = np.linalg.norm(coords[i]-coords[j])
			avgIntermolDist += dist
			nPairs += 1
	avgIntermolDist /= nPairs
	return avgINtermolDist


# from directory name, determine which monomers are in dimer
def getMonNames(file, monNameList): 
	monNameList.sort(key = len, reverse = True)
	name = file.split('/')[1].lower()
	mon1 = 'dummy'
	mon2 = 'dummy'
	for mon in monNameList:
		if mon in name:
			mon1 = mon
			if 'dim' in name:
				mon2 = mon1
			name = name.replace(mon, '')
			break
	if mon2 == 'dummy':
		mon2 = name
	return mon1, mon2

#implement brute-force planarity measure described in Antic2017
# note the the coord array has coords for ALL ATOMS, not just the ones in 
# the monomer that the index set picks out 
def planarity(coordArray): 
	nAtoms = coordArray.shape[0] 
	triples = combinations(range(nAtoms), 3)
	minPlanarity = 100000
	for t in triples: 
		vec1 = coordArray[t[0]] - coordArray[t[1]]
		vec2 = coordArray[t[0]] - coordArray[t[2]] 
		normal = np.cross(vec1, vec2) 
		const = 0 #norm and const define this triple's plane
		for i in range(0,3): 
			const += -1 * normal[i] * coordArray[t[0],i]
		totDist = 0
		for i in range(0,nAtoms): 
			dist = const
			for j in range(0,3): 
				dist += normal[j]*coordArray[i,j]
			dist = abs(dist) / np.linalg.norm(normal) 
			totDist += dist 
		meanDist = totDist / (nAtoms - 3) 
		if meanDist < minPlanarity: 
			minPlanarity = meanDist
	return minPlanarity 

#We define diameter as the maximum C-C distance in the molecule
def diameter(coordArray): 
	nAtoms = coordArray.shape[0]
	pairs = combinations(range(nAtoms), 2) 
	di = 0
	for p in pairs: 
		dist = np.linalg.norm(coords[p[0]]-coords[p[1]])	
		if dist > di: 
			di = dist
	return di

#first process monomers 
outfiles = glob.glob('*Monomer*/*/S1opt/out.complex*') 
atomMassDict = {'H':0, 'C':12.0, 'O':16.0}
monMassDict = {}
monEmissionDict = {}
monOxygenatedDict = {}
monBandgapDict = {}
monBondLenDict = {}
monS1Dict = {}
monNameList = list()
coordFile = open('coords.dat', 'w')
for file in outfiles: 
	name = file.split('/')[1].lower()
	monNameList.append(name) 
	out = open(file, 'r') 
	outLines = out.readlines()
	nAtom = 0
	lineNum = 0
	wavelength = 0
	done = False 
	oxContaining = False
	for index, line in enumerate(outLines): 
		if name not in monBandgapDict: 
			if ('   2.0000   ' in outLines[index] and '   0.0000   ' in outLines[index+1]):
				eHomo = float(outLines[index].split()[2])
				eLumo = float(outLines[index+1].split()[2])
				bandGap = eLumo - eHomo
				bandGap = 10000000/(bandGap*219474.63)
				monBandgapDict[name] = bandGap
		if (name not in monS1Dict) and (done): 
			if 'FINAL SINGLE POINT ENERGY' in line: 
				tokens = line.split()
				monS1Dict[name] = float(tokens[4])
		if 'Number of atoms' in line: 
			tokens = line.split()
			nAtom = int(tokens[4])
		if('FINAL ENERGY' in line): 
			done = True
			lineNum = index+6
			mass = 0
			for j in range(lineNum,lineNum+nAtom): 
				tokens = outLines[j].split()
				mass += atomMassDict[tokens[0]]
				if tokens[0] == 'O': 
					oxContaining = True
			monMassDict[name] = mass
			monOxygenatedDict[name] = oxContaining
			# calculate mass, read in coords of heavy atoms 
			coords = list()
			for j in range(lineNum,lineNum+nAtom):
				 tokens = outLines[j].split()
				 elem = tokens[0]
				 if elem != 'H':
					xyz = [float(i) for i in tokens[1:]]
					coords.append(xyz)
			nHeavyAt = len(coords)
			coords = np.array(coords)
				# a loop over all pairs 
			bondCounter = 0;
			alpha = 98.89
			bondLens = list()
			avgBondLength = 0; 
			for i in range(0, nHeavyAt-1):
				for j in range(i+1, nHeavyAt):
					# overlap 
					 xyDist = np.linalg.norm(coords[i,:2]-coords[j,:2])
					 xyDistCut = 0.1
					 if xyDist < xyDistCut:
						nOverlap += 1
						# monomer membership 
					 bondCutoff = 1.7
					 dist = np.linalg.norm(coords[i]-coords[j])
					 if dist < bondCutoff:
						 avgBondLength += dist;
						 bondCounter += 1;
						 bondLens.append(dist)
			homa = 0
			avgBondLength /= len(bondLens);
			monBondLenDict[name] = avgBondLength
	if done: 
		for index in range(lineNum, len(outLines)-1): 
			if ('   2.0000   ' in outLines[index] and '   0.0000   ' in outLines[index+1]): 
				eHomo = float(outLines[index].split()[2])
				eLumo = float(outLines[index+1].split()[2])
				bandGap = eLumo - eHomo
				#print(str(format(bandGap, '.3f'))+'  ', end = '') 
				# get bandgap in Ha
				bandGap = 10000000/(bandGap*219474.63)
				#monBandgapDict[name] = bandGap
			if 'ELECTRIC DIPOLE MOMENTS' in outLines[index]: 
				tokens = outLines[index+5].split()
				wavelength = tokens[2]	
				monEmissionDict[name] = float(wavelength)


# now process complexes	 
dirList = glob.glob('Batch*[!Monomer]/*/S1opt')
otherDirs = (glob.glob('BatchNewCoro/*/S1opt'))
for d in otherDirs: 
	dirList.append(d)
printHeader()
#otherOutfiles = ['BatchPentaDims/benzogchrysDim/S1opt/out.complexRestartRestart','/panfs/panbox/home/rkrueger/SpectroscopicScan/Noncovalent/BatchPentaDims/dibenzophenDim/OneRootOpt/out.complex']
for dir in dirList: 
	outfiles = glob.glob(dir+'/out.complex*')
	finalEn = 0
	avgIntermolDist = 0; 
	avgBondLength = 0; 
	#if ('BatchPentaDims' in dir): 
	#	outfiles.extend(otherOutfiles) 
	for file in outfiles:
		La = -1
		out = open(file, 'r') 
		outLines = out.readlines()
		# determine whether job finished
		if ('TERMINATED' not in outLines[len(outLines)-2]): 
			continue #if not, skip this file
		# first we figure out which monomers are represented 
		mon1, mon2 = getMonNames(file, monNameList) 
		# if we don't have both monomers, skip this complex 
		if (mon2 not in monNameList) or (mon1 not in monNameList): 
			continue
		print((mon1+'_'+mon2).ljust(25), end='') #col 0
		nAtom = 0
		lineNum = 0
		wavelength = 0
		oxContaining = False
		calcFinished = False 
		diameterSmall = 0
		diameterBig = 0
		for index, line in enumerate(outLines):
			if 'Number of atoms' in line:
				tokens = line.split()
				nAtom = int(tokens[4])
		# everything printed after this is about the final, optimized struct
			if 'STATE  1' in line:
				exciteLine1 = index
			if 'STATE  2' in line:
				exciteLine2 = index
				cutoff = 0.25
				nAbove = 0
				# determine whether La or Lb excitation 
				for q in range(exciteLine1+1, exciteLine2-1):
					tokens = outLines[q].split()
					contribution = float(tokens[4])
					if contribution > cutoff:
						nAbove += 1
				if nAbove > 1:
					La = 0
				else:
					La = 1
			if 'FINAL SINGLE POINT ENERGY' in line:
				tokens = line.split()
				finalEn = float(tokens[4])
			if('FINAL ENERGY' in line):
				calcFinished = True
				lineNum = index+6
				mass = 0
			# calculate mass, read in coords of heavy atoms 
				coords = list()
				elems = list()
				oxContaining = '0'
				for j in range(lineNum,lineNum+nAtom):
					tokens = outLines[j].split()
					elem = tokens[0]
					if elem != 'H': 
						elems.append(elem) 
						xyz = [float(i) for i in tokens[1:]]
						coords.append(xyz)
					mass += atomMassDict[tokens[0]]
					if elem == 'O':
						oxContaining = '1'
				homodimer = '0'
				if mon1 == mon2: 
					homodimer = '1'
				print(homodimer.ljust(3), end='') # col1 
				print(oxContaining.ljust(3), end='') #col 2
				print(str(format(mass, '.2f'))+'  ', end = '') #col 3
				# now examine how many atoms are "overlapping" 
				# also assign atoms to molecules 
				mol1indices = set() 
				mol2indices = set() 
				pairSet = set()
				nOverlap = 0
				mol1indices.add(0) # place the first atom here 
				nHeavyAt = len(coords) 
				coords = np.array(coords) 
				# a loop over all pairs 
				bondCounter = 0; 
				for i in range(0, nHeavyAt-1): 
					for j in range(i+1, nHeavyAt): 
						# overlap 
						xyDist = np.linalg.norm(coords[i,:2]-coords[j,:2])
						xyDistCut = 0.1
						if xyDist < xyDistCut:	 
							nOverlap += 1
						# monomer membership 
						bondCutoff = 1.7
						dist = np.linalg.norm(coords[i]-coords[j])
						if dist < bondCutoff: 
							pairSet.add((i, j)) 
							avgBondLength += dist; 
							bondCounter += 1; 
				avgBondLength /= bondCounter; 
				done = False 
				while (not done): 
					pairAdded = (-1, -1) 
					for pair in pairSet: 
						added = False
						if (pair[0] in mol1indices and pair[1] not in mol1indices):
							mol1indices.add(pair[1])
							added = True
						if (pair[1] in mol1indices and pair[0] not in mol1indices):
							mol1indices.add(pair[0])
							added = True
						if added: 
							pairAdded = pair 
							break 
					if pairAdded[0] ==  -1: #we failed to add one pair tho we tried them all 
						done = True
					else: 
						pairSet.remove(pairAdded) 
				# the leftovers go to mol2 
				for pair in pairSet: 
					if pair[0] not in mol2indices and pair[0] not in mol1indices: 
						mol2indices.add(pair[0])
					if pair[1] not in mol2indices and pair[1] not in mol1indices: 
						mol2indices.add(pair[1])
				print(str(nOverlap).ljust(5), end = '') # col 4
				bigMonIndices = mol1indices
				smallMonIndices = mol2indices
				if len(mol2indices) > len(mol1indices): 
					bigMonIndices = mol2indices
					smallMonIndices	= mol1indices	
				coordsBig = list()
				coordsSmall = list()
				for ind in bigMonIndices: 
					coordsBig.append(coords[ind])
				for ind in smallMonIndices: 
					coordsSmall.append(coords[ind])
				coordsBig = np.array(coordsBig)
				coordsSmall = np.array(coordsSmall) 
				# determine average intermolecular atomic distance
				avgIntermolDist = getIntermolDist(coordsBig, coordsSmall) 
				# monomer planarity in dimer 
				bigPlanarity = planarity(coordsBig) 
				smallPlanarity = planarity(coordsSmall) 
				print(format(bigPlanarity, '.4f')+'  ', end='') # col 5
				print(format(smallPlanarity, '.4f')+'  ', end='') # col 6
				break
		largerMon = mon1 #print some monomer info 
		smallerMon = mon2 
		if monMassDict[mon2] > monMassDict[mon1]: 
			smallerMon = mon1
			largerMon = mon2
		print(format(monMassDict[largerMon], '.2f')+'  ', end = '') # col 7
		print(format(monMassDict[smallerMon], '.2f')+'  ', end = '') # ol 8
		print(format(monBandgapDict[largerMon], '.5f')+'  ', end = '') # col 9
		print(format(monBandgapDict[smallerMon], '.5f')+'  ', end = '') # col 10
		print(format(monEmissionDict[largerMon], '.4f')+'  ', end = '') # col 11
		print(format(monEmissionDict[smallerMon], '.4f')+'  ', end = '') # col 12
		# The dimer fluorescence emission energy 
		wavelength = 0
		for index in range(lineNum, len(outLines)-1):
			if 'ELECTRIC DIPOLE MOMENTS' in outLines[index]:
				tokens = outLines[index+5].split()
				wavelength = float(tokens[2])
				print(str(format(wavelength, '.4f')), end='  ') # col 13
		ha2kj = 2625.5
		# Print remaining stats
		print(format(diameterBig, '.2f')+'  ', end = '')
		print(format(diameterSmall, '.2f')+'  ', end = '')
		print(format(avgIntermolDist, '.4f')+' ', end = '')
		print(format(avgBondLength, '.4f')+' ', end = '')
		print(str(La)+'  ', end = '')
		print(format(monBondLenDict[largerMon], '.5f'), end=' ')
		print(format(monBondLenDict[smallerMon], '.5f'))
		
       
