#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb
#SBATCH --output=output.log
#SBATCH --partition=Orion


python3 master.py



