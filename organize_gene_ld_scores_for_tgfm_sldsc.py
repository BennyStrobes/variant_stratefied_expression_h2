import numpy as np 
import os
import sys
import pdb
import gzip




def get_pseudotissue_names(gtex_pseudotissue_file):
	f = open(gtex_pseudotissue_file)
	arr = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(data[0])
	f.close()
	return np.asarray(arr)

def get_n_genes_per_chromosome_per_tissue(pseudotissue_names, gtex_susie_gene_models_dir):
	n_genes = np.zeros((22, len(pseudotissue_names)))
	for tissue_iter, pseudotissue_name in enumerate(pseudotissue_names):
		pos_file = gtex_susie_gene_models_dir + pseudotissue_name + '/' + pseudotissue_name + '_cis_heritable_gene_pos_file.txt'
		tmp = np.loadtxt(pos_file, dtype=str,delimiter='\t')
		for chrom_num in range(1,23):
			gene_count = np.sum(tmp[1:,2] == str(chrom_num))
			n_genes[(chrom_num-1), tissue_iter] = gene_count
	return n_genes

def load_in_gene_ld_scores(gene_ld_score_file):
	gene_ld_scores = np.loadtxt(gene_ld_score_file)
	return gene_ld_scores

def merge_variant_and_gene_ld_score_files(variant_ld_score_file, gene_ld_score_file, pseudotissue_name, merged_ld_score_file):
	gene_ld_scores = load_in_gene_ld_scores(gene_ld_score_file)
	gene_ld_scores[np.isnan(gene_ld_scores)] = 0.0

	f = gzip.open(variant_ld_score_file)
	t = open(merged_ld_score_file,'w')
	head_count = 0
	line_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data)))
			for ele in np.asarray(data[3:]):
				t.write('\t' + pseudotissue_name + '_eQTL_' + ele)
			t.write('\n')
			continue
		t.write('\t'.join(np.asarray(data)) + '\t' + '\t'.join(gene_ld_scores[line_counter,:].astype(str)) + '\n')
		line_counter = line_counter + 1
	f.close()
	t.close()

	# Quick error checks before quiting
	if line_counter != gene_ld_scores.shape[0]:
		print('fundamental assumption eroror')
		pdb.set_trace()
	return 

def merge_variant_and_gene_ld_score_files_for_genotype_intercept(variant_ld_score_file, gene_ld_score_file, pseudotissue_name, merged_ld_score_file):
	gene_ld_scores = load_in_gene_ld_scores(gene_ld_score_file)
	gene_ld_scores[np.isnan(gene_ld_scores)] = 0.0

	f = gzip.open(variant_ld_score_file)
	t = open(merged_ld_score_file,'w')
	head_count = 0
	line_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(np.asarray(data)) + '\t' + pseudotissue_name + '_eQTL_baseL2' + '\n')
			continue
		t.write('\t'.join(np.asarray(data)) + '\t' + str(gene_ld_scores[line_counter,0]) + '\n')
		line_counter = line_counter + 1
	f.close()
	t.close()

	# Quick error checks before quiting
	if line_counter != gene_ld_scores.shape[0]:
		print('fundamental assumption eroror')
		pdb.set_trace()
	return 

def merge_m_files(variant_m_file, gene_m_file, merged_m_file):
	gene_m_vec = np.loadtxt(gene_m_file)
	t = open(merged_m_file,'w')
	m_vec = np.loadtxt(variant_m_file)
	t.write('\t'.join(m_vec.astype(str)) + '\t' + '\t'.join(gene_m_vec.astype(str)) + '\n')
	t.close()

def merge_m_files_for_genotype_intercept(variant_m_file, gene_m_file, merged_m_file):
	gene_m_vec = np.loadtxt(gene_m_file)

	t = open(merged_m_file,'w')
	m_vec = np.loadtxt(variant_m_file)
	t.write('\t'.join(m_vec.astype(str)) + '\t' + str(gene_m_vec[0]) + '\n')
	t.close()

def merge_variant_and_gene_annotation_files(variant_annot_file, gene_annotation_file, output_anno_file, pseudotissue_name, chrom_num):
	t = open(output_anno_file,'w')
	f = open(variant_annot_file)
	head_count = 0
	# Stream variant annotation file
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[4:])
			header = np.asarray(data[:4])
			t.write('\t'.join(header) + '\t' + '\t'.join(anno_names))
			for anno_name in anno_names:
				t.write('\t' + pseudotissue_name + '_eQTL_' + anno_name)
			t.write('\n')
			num_variant_annotations = len(anno_names)
			num_eqtl_annotations = len(anno_names)
			continue
		t.write('\t'.join(np.asarray(data)) + '\t' + '\t'.join(np.zeros(num_eqtl_annotations).astype(str)) + '\n')
	f.close()
	# Stream gene annotation file
	f = open(gene_annotation_file)
	placeholder_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Print header columns
		t.write(str(chrom_num) + '\t' + str(placeholder_counter) + '\t' + data[0] + '\t' + str(placeholder_counter))
		t.write('\t' + '\t'.join(np.zeros(num_variant_annotations).astype(str)))
		t.write('\t' + '\t'.join(np.asarray(data[1:])) + '\n')
		placeholder_counter = placeholder_counter + 1
	f.close()
	t.close()

def merge_variant_and_gene_intercept_annotation_files(variant_annot_file, gene_annotation_file, output_anno_file, pseudotissue_name, chrom_num):
	t = open(output_anno_file,'w')
	f = open(variant_annot_file)
	head_count = 0
	# Stream variant annotation file
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[4:])
			header = np.asarray(data[:4])
			t.write('\t'.join(header) + '\t' + '\t'.join(anno_names))
			for anno_name in anno_names[0:1]:
				t.write('\t' + pseudotissue_name + '_eQTL_' + anno_name)
			t.write('\n')
			num_variant_annotations = len(anno_names)
			num_eqtl_annotations = 1
			continue
		t.write('\t'.join(np.asarray(data)) + '\t' + '\t'.join(np.zeros(num_eqtl_annotations).astype(str)) + '\n')
	f.close()
	# Stream gene annotation file
	f = open(gene_annotation_file)
	placeholder_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Print header columns
		t.write(str(chrom_num) + '\t' + str(placeholder_counter) + '\t' + data[0] + '\t' + str(placeholder_counter))
		t.write('\t' + '\t'.join(np.zeros(num_variant_annotations).astype(str)))
		t.write('\t' + '\t'.join(np.asarray(data[1:2])) + '\n')
		placeholder_counter = placeholder_counter + 1
	f.close()
	t.close()


def merge_variant_and_gene_frq_files(input_frq_file, gene_annotation_file, output_frq_file, chrom_num):
	t = open(output_frq_file,'w')
	f = open(input_frq_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		t.write('\t'.join(np.asarray(data)) + '\n')
	f.close()
	f = open(gene_annotation_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		t.write(str(chrom_num) + '\t' + data[0] + '\tA\tG\t0.25\t978\n')
	f.close()
	t.close()



def make_ld_score_input_shell(chrom_num, pseudotissue_name, variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, ref_1kg_genotype_dir):
	input_frq_file = ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.frq'
	for variant_model in variant_models:
		variant_annot_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.annot'
		variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
		variant_m_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.M'
		variant_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + variant_model + '.' + str(chrom_num) + '.l2.M_5_50'
		for gene_model_suffix in gene_model_suffixes:
			gene_annotation_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + 'gene_annotations.txt'
			gene_ld_score_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + gene_model_suffix
			gene_m_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_gene_m'
			gene_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_gene_m_5_50'

			# Create new frq file
			output_frq_file = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_' + variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.frq'
			merge_variant_and_gene_frq_files(input_frq_file, gene_annotation_file, output_frq_file, chrom_num)


			# Create new annotation file
			output_anno_file = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_' + variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.annot'
			merge_variant_and_gene_annotation_files(variant_annot_file, gene_annotation_file, output_anno_file, pseudotissue_name, chrom_num)

			# Create ld score output root
			output_root = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_' + variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.l2'
			print(variant_model + '  ' + gene_model_suffix)
			# Merge ld scores files
			merged_ld_score_file = output_root + '.ldscore'
			merge_variant_and_gene_ld_score_files(variant_ld_score_file, gene_ld_score_file, pseudotissue_name, merged_ld_score_file)

			merged_m_file = output_root + '.M'
			merge_m_files(variant_m_file, gene_m_file, merged_m_file)

			merged_m_5_50_file = output_root + '.M_5_50'
			merge_m_files(variant_m_5_50_file, gene_m_5_50_file, merged_m_5_50_file)

def make_ld_score_input_shell_for_genotype_intercept(chrom_num, pseudotissue_name, reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir):
	input_frq_file = ref_1kg_genotype_dir + '1000G.EUR.hg38.' + str(chrom_num) + '.frq'
	variant_ld_score_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.ldscore.gz'
	variant_m_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.M'
	variant_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.l2.M_5_50'
	variant_annot_file = preprocessed_tgfm_sldsc_data_dir + reference_variant_model + '.' + str(chrom_num) + '.annot'
	for gene_model_suffix in gene_model_suffixes:
		gene_annotation_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + 'gene_annotations.txt'
		gene_ld_score_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_' + gene_model_suffix
		gene_m_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_gene_m'
		gene_m_5_50_file = preprocessed_tgfm_sldsc_data_dir + 'tissue_eqtl.' + str(chrom_num) + '.' + pseudotissue_name + '_gene_m_5_50'

		# Create new frq file
		output_frq_file = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_eqtl_intercept_' + reference_variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.frq'
		merge_variant_and_gene_frq_files(input_frq_file, gene_annotation_file, output_frq_file, chrom_num)
		
		# Create new annotation file
		output_anno_file = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_eqtl_intercept_' + reference_variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.annot'
		merge_variant_and_gene_intercept_annotation_files(variant_annot_file, gene_annotation_file, output_anno_file, pseudotissue_name, chrom_num)

		# Create output root
		output_root = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_eqtl_intercept_' + reference_variant_model + '_' + gene_model_suffix + '.' + str(chrom_num) + '.l2'

		print('eqtl_intercept' + '  ' + gene_model_suffix)
		# Merge ld scores files
		merged_ld_score_file = output_root + '.ldscore'
		merge_variant_and_gene_ld_score_files_for_genotype_intercept(variant_ld_score_file, gene_ld_score_file, pseudotissue_name, merged_ld_score_file)

		merged_m_file = output_root + '.M'
		merge_m_files_for_genotype_intercept(variant_m_file, gene_m_file, merged_m_file)

		merged_m_5_50_file = output_root + '.M_5_50'
		merge_m_files_for_genotype_intercept(variant_m_5_50_file, gene_m_5_50_file, merged_m_5_50_file)

def save_genetic_element_partition(output_genetic_element_parition_file, partition_names, partition_positions):
	t = open(output_genetic_element_parition_file, 'w')
	t.write('partition_name\tpartition_start\tpartition_end\n')
	for indexer,partition_name in enumerate(partition_names):
		t.write(partition_name + '\t' + str(partition_positions[indexer][0]) + '\t' + str(partition_positions[indexer][1]) + '\n')
	t.close()


gtex_pseudotissue_file = sys.argv[1]
preprocessed_tgfm_sldsc_data_dir = sys.argv[2]
gene_type = sys.argv[3]
gtex_susie_gene_models_dir = sys.argv[4]
ref_1kg_genotype_dir = sys.argv[5]


# Extract names of pseudotissues
#pseudotissue_names = get_pseudotissue_names(gtex_pseudotissue_file)
pseudotissue_name = 'Whole_Blood'

# Get number of genes per chromosome per tissue
#n_genes_per_chromosome_per_tissue = get_n_genes_per_chromosome_per_tissue(pseudotissue_names, gtex_susie_gene_models_dir)

# Various iterations to run over
variant_models = ['baselineLD_no_qtl']

# Gene modedls

gene_model_suffixes = ['gene_ld_scores', 'gene_adj_ld_scores']
for chrom_num in range(1,23):
	print(chrom_num)
	make_ld_score_input_shell(chrom_num,  pseudotissue_name, variant_models, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir, ref_1kg_genotype_dir)

# Make version without genomic annotaitons (call it genotype intercept)
reference_variant_model = 'baselineLD_no_qtl'
for chrom_num in range(1,23):
	make_ld_score_input_shell_for_genotype_intercept(chrom_num, pseudotissue_name, reference_variant_model, gene_model_suffixes, preprocessed_tgfm_sldsc_data_dir)



# Save variant partitions
for variant_model in variant_models:
	for gene_model_suffix in gene_model_suffixes:
		output_genetic_element_parition_file = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_' + variant_model + '_' + gene_model_suffix + '.genetic_element_partition.txt'
		save_genetic_element_partition(output_genetic_element_parition_file, ['non-mediated-variant', pseudotissue_name + '-eqtl-component'],[[0,92], [93, 185]])

# No genomic annotations
for gene_model_suffix in gene_model_suffixes:
	output_genetic_element_parition_file = preprocessed_tgfm_sldsc_data_dir + pseudotissue_name + '_eqtl_intercept_' + reference_variant_model + '_' + gene_model_suffix + '.genetic_element_partition.txt'
	save_genetic_element_partition(output_genetic_element_parition_file, ['non-mediated-variant', pseudotissue_name + '-eqtl-component'],[[0,92], [93, 93]])



