#!/bin/bash
#SBATCH --job-name=extract_timepoints
#SBATCH --output=/home/lshpk18/extract_timepoints_%j.log
#SBATCH --error=/home/lshpk18/extract_timepoints_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate dualrnaseq_deseq2

bash "/home/lshpk18/unmapped reads (kraken)/extract_by_timepoint.sh"
