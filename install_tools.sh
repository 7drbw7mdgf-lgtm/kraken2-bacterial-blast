#!/bin/bash
set -e

echo "Installing Kraken2 into base environment..."
conda install -y -c bioconda -c conda-forge kraken2
echo "Kraken2 installed: $(kraken2 --version | head -1)"

echo "Creating unicycler environment with Python 3.12 (base env has Python 3.13 which is incompatible)..."
conda create -y -n unicycler -c bioconda -c conda-forge python=3.12 unicycler
echo "Unicycler installed in 'unicycler' conda environment."
echo "To use Unicycler: conda activate unicycler"

echo "All tools installed successfully."
