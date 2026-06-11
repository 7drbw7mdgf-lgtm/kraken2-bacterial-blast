#!/bin/bash
#SBATCH --job-name=extract_cutibacterium
#SBATCH --output=/home/lshpk18/cutibacterium_extract_%j.log
#SBATCH --error=/home/lshpk18/cutibacterium_extract_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

bash "/home/lshpk18/unmapped reads (kraken)/extract_cutibacterium_reads.sh"
