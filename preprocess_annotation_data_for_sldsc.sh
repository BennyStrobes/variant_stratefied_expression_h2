#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20GB                         # Memory total in MiB (for all cores)




ldsc_code_dir="$1"
hapmap3_rsid_file="$2"
ldsc_baseline_annotation_dir="$3"
ldsc_baseline_ld_annotation_dir="$4"
ref_1kg_genotype_dir="$5"
gtex_pseudotissue_file="$6"
gtex_susie_gene_models_dir="$7"
preprocessed_tgfm_sldsc_data_dir="$8"
gene_type="${9}"
chrom_num="${10}"


echo $chrom_num

##################################
# Create variant annotation files
###################################
# input file (baselineld annotation file)
baseline_ld_annotation_stem=${ldsc_baseline_ld_annotation_dir}baselineLD.${chrom_num}
# output file (baselineld annotation file without any qtl maxcpp annotations)
baseline_ld_no_qtl_annotation_stem=${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}
# Perform filtering
source ~/.bash_profile
python3 remove_qtl_annotations_from_baselineld_annotation_file.py $baseline_ld_annotation_stem $baseline_ld_no_qtl_annotation_stem


# Create baseline annotation file
cp ${ldsc_baseline_annotation_dir}baseline.${chrom_num}.annot.gz ${preprocessed_tgfm_sldsc_data_dir}baseline_no_qtl.${chrom_num}.annot.gz



##################################
# Create variant level LD-scores
###################################
# Load in LDSC module
source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12

# Run standard S-LDSC on baselineLD annotations
python ${ldsc_code_dir}ldsc.py\
	--l2\
	--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
	--ld-wind-cm 1\
	--annot ${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}.annot\
	--out ${preprocessed_tgfm_sldsc_data_dir}baselineLD_no_qtl.${chrom_num}\
	--print-snps ${hapmap3_rsid_file}
# Run standard S-LDSC on baseline annotations
python ${ldsc_code_dir}ldsc.py\
	--l2\
	--bfile ${ref_1kg_genotype_dir}1000G.EUR.hg38.${chrom_num}\
	--ld-wind-cm 1\
	--annot ${preprocessed_tgfm_sldsc_data_dir}baseline_no_qtl.${chrom_num}.annot.gz\
	--out ${preprocessed_tgfm_sldsc_data_dir}baseline_no_qtl.${chrom_num}\
	--print-snps ${hapmap3_rsid_file}




