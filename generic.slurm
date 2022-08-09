#!/bin/bash

#SBATCH --job-name=ASE
#SBATCH --mincpus=8
#SBATCH --mem-per-cpu=6GB
#SBATCH --time=4:00:00
#SBATCH --partition=shared

cd /gpfs/home/elu/ASEClass/

module load python/3.8.3

# python3 transcriptFilter.py
python3 ExonClassifier.py
