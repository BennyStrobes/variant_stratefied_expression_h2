#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40G                         # Memory total in MiB (for all cores)




sldsc_results_dir="$1"


source ~/.bash_profile


python3 meta_analyze_sldsc_results_for_blood_cell_traits.py $sldsc_results_dir