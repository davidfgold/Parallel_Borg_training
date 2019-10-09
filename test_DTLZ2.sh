#!/bin/bash
#SBATCH --job-name="dtlz2_1"
#SBATCH --output="dtlz2_01.out"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --export=ALL
#SBATCH -t 0:1:00

#ibrun in verbose mode will give binding detail
#this runs UF11 test problem for 1 minute

mpirun -n 128 ./dtlz2_ms.exe 1