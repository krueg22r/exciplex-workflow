#! /usr/bin/env python3

# set of functions to write scripts for Maui/Torque queueing system on 
# an HPC cluster for quantum chemistry packages Orca and Gaussian 

import math
import re

#file path info
mpiPath = '/panfs/panbox/home/krueger/openmpi2/'
orcaExe = '/panfs/panbox/home/rkrueger/orca-4.0/orca'

#prints info for scheduler, accepts open file for writing
def printHeader(sh, nProc): 
	print("#!/bin/bash", file = sh)
	print("#PBS -V",  file = sh)
	print("#PBS -M rkrueger@caltech.edu", file = sh)
        print("#PBS -m ae", file = sh)
        nNode = math.ceil(int(nProc) / 12)
	if(int(nProc) >= 12):
                print("#PBS -l nodes="+str(nNode)+":ppn=12", file = sh)
        if(int(nProc) < 12): print("#PBS -l nodes="+str(nNode)+":ppn="+str(nProc), file = sh)

#prints info to update environment variables
def printEnvInfo(sh): 
	print("export RSH_COMMAND=\"/usr/bin/ssh -x\" ", file = sh)
        print("export PATH="+mpiPath+"bin:$PATH", file=sh)
        print("export LD_LIBRARY_PATH="+mpiPath+"lib:$LD_LIBRARY_PATH",file=sh)
	print('export TMPDIR=/tmp', file=sh)
        print('cd $PBS_O_WORKDIR', file=sh)

# For Orca QC code, 330 hour limit 
def writeScriptLong(baseName, nProc=12):
        shName = baseName+".sh"
        sh=open(shName, 'w')
	printHeader(sh, nProc)
        print("#PBS -N "+baseName, file = sh)
        print("#PBS -l walltime=330:00:00", file=sh)
        print("#PBS -q fast", file = sh)
        print("", file = sh)
	printEnvInfo(sh) 
        print(orcaExe+" $PBS_O_WORKDIR/"+baseName+" > $PBS_O_WORKDIR/out."+baseName, file = sh)

# For Orca QC code, 48-hour limit 
def writeScriptShort(baseName, nProc=12):
        shName = baseName+".sh"
        sh=open(shName, 'w')
	printHeader(sh, nProc) 
        print("#PBS -N "+baseName, file = sh)
        print("#PBS -l walltime=48:00:00", file=sh)
        print("#PBS -q default", file = sh)
        print("", file = sh)
	printEnvInfo(sh) 
        print("$orcadir/orca $PBS_O_WORKDIR/"+baseName+" > $PBS_O_WORKDIR/out."+baseName, file = sh)

# for another quantum chemistry code, Gaussian 
def writeGaussianScript(baseName, nProc): 
        shName = baseName+".sh"
        sh=open(shName, 'w')
	printHeader(sh, nProc) 
        print("#PBS -l walltime=48:00:00", file = sh)
        print("#PBS -q default", file = sh)
        print("", file = sh)
        print("cd $PBS_O_WORKDIR", file = sh)
        print("G_DIR=/panfs/panbox/home/rkrueger/g09", file = sh)
        print("$G_DIR/g09 < "+baseName+" > out."+baseName, file = sh)
