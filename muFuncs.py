#! /usr/bin/env python3

# this function calculates the value of the range-separation 
# parameter mu that minimizes the error function J^2, which 
# derives from the DFT version of Koopman's theorem 

import glob
import numpy as np 
from scipy import interpolate 

# read Orca output file, extract energy info + the mu used for this run
def procFile(fileName):
        f = open(fileName, 'r')
        out = f.readlines()
        eHomo = 0
        eTot = 0
        mu = 0
        for i, line in enumerate(out):
                
                if 'FINAL SINGLE POINT ENERGY' in line:
                        tokens = line.split()
                        eTot = float(tokens[4])
                if ('   2.0000   ' in line or '   1.0000   ' in line) and '   0.0000   ' in out[i+1]: 
                        tokens = line.split()
                        eHomo = float(tokens[2])
                if 'RangeSepMu' in line: 
                        tokens = line.split()
                        mu = float(tokens[3])
        return eTot, eHomo, mu

def getMu(dirName): 
        dirList = glob.glob(dirName+'/0.*')
        j2list = list()
        mulist = list()
        for d in dirList:
                nFile = d+'/0/out.mon'
                pFile = d+'/1/out.mon'
		# process neutral mol output
                nTot, nHomo, mu = procFile(nFile)
		#process positive mol output
                pTot, pHomo, mu = procFile(pFile)
                nIP = pTot - nTot
		#calculate error
                j2 = (nHomo+nIP)*(nHomo+nIP)
                j2list.append(j2)
                mulist.append(mu) 
	# use spline fit to calculate mu corresponding to lowest error	
        j2list = np.asarray(j2list)
        mulist = np.asarray(mulist) 
        mulist, j2list = zip(*sorted(zip(mulist, j2list)))
        xVals = np.linspace(np.amin(mulist)-0.1, np.amax(mulist)+0.1, 220) 
        tck = interpolate.splrep(mulist, j2list, s=0, k=2)
        yVals = interpolate.splev(xVals, tck, der=0)
        mu = xVals[np.argmin(yVals)]
        return mu

