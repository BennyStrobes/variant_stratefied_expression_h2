##################
# Input data
##################

# Ldscore regression code
ldsc_code_dir="/n/groups/price/ben/code_temp4/modified_sldsc/"

# Hapmap3 rsids
hapmap3_rsid_file="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/w_hm3.noMHC.snplist"

# LDSC baseline LD Dir
ldsc_baseline_ld_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baselineLD_v2.2/"
ldsc_baseline_annotation_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/baseline_v1.2/"

# Genotype data from 1KG
ref_1kg_genotype_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

# File containing gtex tissues to do analysis on and their sample size
gtex_pseudotissue_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_meta_analysis_eqtl_calling/pseudotissue_sample_names/pseudotissue_info.txt"

# Directory containing GTEx Susie gene models
gtex_susie_gene_models_dir="/n/scratch3/users/b/bes710/causal_eqtl_gwas/gtex/gtex_susie_gene_models/"


# Summary statistics
full_sumstat_dir="/n/groups/price/ldsc/sumstats_formatted_2021/"

# hg38 sldsc weights
sldsc_h38_weights_dir="/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3_hg38/weights/"

##################
# Output data
##################
# Root node
output_root="/n/scratch3/users/b/bes710/variant_stratefied_expression_h2/"

# Directory containing preprocessed SLDSC data
preprocessed_sldsc_annotation_data_dir=$output_root"preprocessed_sldsc_annotation_data/"


sldsc_results_dir=$output_root"sldsc_results/"




########################################
# Preprocess annotation data for S-LDSC
########################################
gene_type="cis_heritable_gene"
if false; then
for chrom_num in $(seq 1 22); do 
	sbatch preprocess_annotation_data_for_sldsc.sh $ldsc_code_dir $hapmap3_rsid_file $ldsc_baseline_annotation_dir $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $gtex_pseudotissue_file $gtex_susie_gene_models_dir $preprocessed_sldsc_annotation_data_dir $gene_type $chrom_num
done
fi

if false; then
sed 1d $gtex_pseudotissue_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	sbatch preprocess_gene_ld_scores_for_tgfm_sldsc.sh $ldsc_code_dir $hapmap3_rsid_file $ldsc_baseline_annotation_dir $ldsc_baseline_ld_annotation_dir $ref_1kg_genotype_dir $pseudotissue_name $gtex_susie_gene_models_dir $preprocessed_sldsc_annotation_data_dir $gene_type
done
fi
if false; then
sh organize_gene_ld_scores_for_tgfm_sldsc.sh $gtex_pseudotissue_file $preprocessed_sldsc_annotation_data_dir $gene_type $gtex_susie_gene_models_dir $ref_1kg_genotype_dir
fi

trait_name="blood_RBC_DISTRIB_WIDTH"
if false; then
sbatch run_tgfm_sldsc.sh $preprocessed_sldsc_annotation_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $sldsc_results_dir $trait_name

trait_name="blood_RED_COUNT"
sbatch run_tgfm_sldsc.sh $preprocessed_sldsc_annotation_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $sldsc_results_dir $trait_name

trait_name="blood_WHITE_COUNT"
sbatch run_tgfm_sldsc.sh $preprocessed_sldsc_annotation_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $sldsc_results_dir $trait_name

trait_name="blood_PLATELET_COUNT"
sbatch run_tgfm_sldsc.sh $preprocessed_sldsc_annotation_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $sldsc_results_dir $trait_name

trait_name="blood_EOSINOPHIL_COUNT"
sbatch run_tgfm_sldsc.sh $preprocessed_sldsc_annotation_data_dir $full_sumstat_dir $ldsc_code_dir $sldsc_h38_weights_dir $ref_1kg_genotype_dir $sldsc_results_dir $trait_name
fi


sh meta_analyze_sldsc_results_for_blood_cell_traits.sh $sldsc_results_dir









