#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)








gtex_pseudotissue_file="$1"
preprocessed_tgfm_sldsc_data_dir="$2"
gene_type="$3"
gtex_susie_gene_models_dir="$4"
ref_1kg_genotype_dir="$5"

source ~/.bash_profile


python3 organize_gene_ld_scores_for_tgfm_sldsc.py $gtex_pseudotissue_file $preprocessed_tgfm_sldsc_data_dir $gene_type $gtex_susie_gene_models_dir $ref_1kg_genotype_dir