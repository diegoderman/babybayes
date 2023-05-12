#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=dderman@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=0-4:59:00
#SBATCH --mem=256gb
#SBATCH --partition=gpu
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=dvars
#SBATCH --output=dvars_%j.out
#SBATCH --error=dvars_%j.err

######  Module commands #####


######  Job commands go below this line #####

python3 ./1.1_calculate-dvars.py --in ../data/table.csv --out ../data/table_dvars.csv --threads 48
python3 ./1.2_filter-dvars.py 
