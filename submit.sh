#!/bin/bash
#SBATCH -J gpu_test
#SBATCH -n 2
#SBATCH -N 2
#SBATCH -c 4
#SBATCH -w cn15,cn16
#SBATCH -t 00:30:00
#SBATCH -o out
#SBATCH -p gpu
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact ./a.out 1200 1200 300 1
