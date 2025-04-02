#!/bin/bash
#SBATCH -p main
#SBATCH -n4
#SBATCH --time=10:00:00

source /scratch/lustre/home/aust0075/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
export PATH=/scratch/lustrte/home/aust0075/miniforge3/bin:$PATH
snakemake -p --use-conda --cores 4
