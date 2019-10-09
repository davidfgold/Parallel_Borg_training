#!/bin/bash
#SBATCH --job-name="UF11_benchmark_9"
#SBATCH --output="UF11_01.out"
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --export=ALL
#SBATCH -t 0:3:00

mpirun -n 128 ./uf11_mm.exe 1