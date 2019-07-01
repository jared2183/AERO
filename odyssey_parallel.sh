#!/bin/bash
#SBATCH -n 24               # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-00:10          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test  	    # Partition to submit to
#SBATCH --mem=5000          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid

module load Anaconda3/5.0.1-fasrc02  #Load Anaconda module
source activate path_planning
time python fleury_general.py 24
