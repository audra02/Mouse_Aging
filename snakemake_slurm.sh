#!/bin/bash
#SBATCH --job-name=chipseq_pipeline
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --output=logs/slurm/%x_%j.out
#SBATCH --error=logs/slurm/%x_%j.err
#SBATCH --mem=8G

# Load necessary modules and conda
source /scratch/lustre/home/mekl0332/miniforge3/etc/profile.d/conda.sh
conda activate snakemake_env
export PATH=/scratch/lustre/home/mekl0332/miniforge3/bin:$PATH

# Create log directory if it doesn't exist
mkdir -p logs/slurm

# Run Snakemake with SLURM profile
snakemake -p \
    --use-conda \
    --conda-frontend conda \
    --cores $SLURM_CPUS_PER_TASK \
    --latency-wait 60 \
    --restart-times 1 \
    --keep-going \
    --jobs 4
