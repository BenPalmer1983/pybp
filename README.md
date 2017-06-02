# pybp


To run in the command line:

Set the number of processes (otherwise, a default of 4 will be used):


user  /home $      python3 /pathtofile/converge.py example.in



To run with MOAB:

#!/bin/bash
#########################################################
#
# An example script to run the script on BlueBEAR
# Loads the required python modules
# Loads QE 6.0 module
# Sets procCount environment variable (although this is 
# set again in the python script).
# Runs the python script
#
# PWscf batch file
# Settings for Moab job scheduler
#MOAB -l "nodes=4:ppn=8,walltime=4:00:00"
#MOAB -j oe
#MOAB -A readmsd02
#
# Load python3
module load apps/python3/v3.3.0
module load apps/openblas/v0.2.11_gcc-4.7.2
# Load the QE6.0 module
module load espresso/gcc/6.0
# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"
# Set the number of threads to 1
export OMP_NUM_THREADS=1

python3 $PYAPP/converge.py example.in   



The modules are specific to the computer that I wrote this batch script for, and $PYAPP is an environment variable I set on that computer, and it stores the path to my python scripts.







