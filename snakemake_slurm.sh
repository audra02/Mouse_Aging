#!/bin/bash
#SBATCH --job-name=chipseq
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output=logs/slurm/%x_%j.out
#SBATCH --error=logs/slurm/%x_%j.err

# Clean problematic directories
rm -rf results/trimmed/ results/aligned/ logs/trim_galore/

source /scratch/lustre/home/mekl0332/miniforge3/etc/profile.d/conda.sh
conda activate snakemake_env
export PATH=/scratch/lustre/home/mekl0332/miniforge3/bin:$PATH

# Create log directory if it doesn't exist
mkdir -p logs/slurm

export PATH=/scratch/lustre/home/aust0075/miniforge3/bin:$PATH
# Main execution
snakemake --unlock
snakemake -p \
    --use-conda \
    --conda-frontend conda \
    --cores $SLURM_CPUS_PER_TASK \
    --latency-wait 60 \
    --restart-times 1 \
    --keep-going \
    --jobs 4
