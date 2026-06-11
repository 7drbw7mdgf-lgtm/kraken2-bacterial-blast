#!/bin/bash
#SBATCH --job-name=batch_kraken2
#SBATCH --output=/home/lshpk18/kraken2_slurm_%j.log
#SBATCH --error=/home/lshpk18/kraken2_slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

bash "/home/lshpk18/unmapped reads/batch_kraken2.sh"
