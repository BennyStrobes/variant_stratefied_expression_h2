#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40G                         # Memory total in MiB (for all cores)




preprocessed_tgfm_sldsc_data_dir="$1"
full_sumstat_dir="$2"
ldsc_code_dir="$3"
sldsc_h38_weights_dir="$4"
ref_1kg_genotype_dir="$5"
tgfm_sldsc_results_dir="$6"
trait_name="$7"


source /n/groups/price/ben/environments/sldsc/bin/activate
module load python/2.7.12


trait_file=$full_sumstat_dir"UKB_460K."$trait_name".sumstats"
data_version="Whole_Blood_baselineLD_no_qtl_gene_adj_ld_scores"
if false; then
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --genetic-element-partitions ${preprocessed_tgfm_sldsc_data_dir}${data_version}".genetic_element_partition.txt" --frqfile-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
fi
python organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}"."



trait_file=$full_sumstat_dir"UKB_460K."$trait_name".sumstats"
data_version="Whole_Blood_eqtl_intercept_baselineLD_no_qtl_gene_adj_ld_scores"
if false; then
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --genetic-element-partitions ${preprocessed_tgfm_sldsc_data_dir}${data_version}".genetic_element_partition.txt" --frqfile-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
fi
python organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}"."














if false; then
# OLD W/O anno-hack
data_version="Whole_Blood_baselineLD_no_qtl_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
python organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}"."

data_version="Whole_Blood_eqtl_intercept_baselineLD_no_qtl_gene_adj_ld_scores"
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${preprocessed_tgfm_sldsc_data_dir}${data_version}"." --w-ld-chr ${sldsc_h38_weights_dir}"weights.hm3_noMHC." --print-coefficients --out ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_"
python organize_tgfm_sldsc_results.py ${tgfm_sldsc_results_dir}${trait_name}"_"${data_version}"_" ${preprocessed_tgfm_sldsc_data_dir}${data_version}"."
fi







if false; then
/n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files
python ${ldsc_code_dir}ldsc.py --h2 ${trait_file} --ref-ld-chr ${ldsc_baseline_ld_hg19_annotation_dir}"baselineLD." --w-ld-chr ${ldsc_weights_dir}"weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${ldsc_genotype_dir}"1000G.EUR.QC." --out ${standard_sldsc_results}${trait_name}"_sldsc_res_"
fi