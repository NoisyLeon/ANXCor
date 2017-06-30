#!/bin/bash

#SBATCH -J S2Cpp
#SBATCH -o S2Cpp_%j.out
#SBATCH -e S2Cpp_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --mem=MaxMemPerNode

module load gcc/gcc-4.9.1 
module load fftw
#module load fftw/fftw-3.3.3_openmpi-1.7.4_gcc-4.8.2_double_ib
#. ~/.my.bashforSEED2COR


cppexe=/projects/life9360/code/Programs/SEED2CORpp_1.08/Seed2Cor
#cd /lustre/janus_scratch/life9360/COR_TEST_forye_Juan_not_syncnorm
cd /lustre/janus_scratch/life9360/COR_TEST_forye_Juan_002
#cd /lustre/janus_scratch/life9360/COR_US_PO_new
#cd /lustre/janus_scratch/life9360/TEST_001
export OMP_NUM_THREADS=12; $cppexe parameters.txt <<- END
Y
END
