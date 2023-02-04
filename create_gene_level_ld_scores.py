import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import os
import pdb
import numpy as np
from pandas_plink import read_plink1_bin
import pickle
import pandas as pd
import pyreadr
import gzip
import time




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



def create_tissue_to_gene_model_df(pseudotissue_name, gtex_susie_gene_models_dir, gene_type, chrom_num):
	gene_model_summary_file = gtex_susie_gene_models_dir + pseudotissue_name + '/' + pseudotissue_name + '_' + gene_type + '_pos_file.txt'
	raw_table = np.loadtxt(gene_model_summary_file, dtype=str, delimiter='\t')[1:]
	rows_on_chromosome = raw_table[:,2] == chrom_num
	raw_table_chrom = raw_table[rows_on_chromosome, :]
	return raw_table_chrom



def create_mapping_from_rsid_to_snpid(rsids, snp_ids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		if rsid in mapping:
			print('assumption eroror')
			pdb.set_trace()
		mapping[rsid] = snp_ids[ii]
	return mapping

def create_mapping_from_snpid_to_rsid(rsids, snp_ids):
	mapping = {}
	for ii, rsid in enumerate(rsids):
		if snp_ids[ii] in mapping:
			print('assumption eororor')
			pdb.set_trace()
		mapping[snp_ids[ii]] = rsid
	return mapping

def create_mapping_snpid_to_reference_index(snp_ids, alt_snp_ids):
	mapping = {}
	for ii, snpid in enumerate(snp_ids):
		mapping[snpid] = (ii, 1.0)
		alt_snpid = alt_snp_ids[ii]
		if alt_snpid == snpid:
			print('assumption eroror')
			pdb.set_trace()
		mapping[alt_snpid] = (ii, -1.0)

	return mapping

def get_gene_model_snp_positions(variant_names):
	positions = []
	for variant_name in variant_names:
		positions.append(variant_name.split('_')[1])
	return np.asarray(positions).astype(int)

def project_gene_model_onto_full_window(gene_susie_mu, gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=True):
	# First get number of components
	n_components = gene_susie_mu.shape[0]

	# Initialize gene model on full window
	window_susie_mu = np.zeros((n_components, n_window_snps))

	# Loop through gene snpid names
	for ii, gene_snpid_name in enumerate(gene_snpid_names):
		# Get window position and sign
		if gene_snpid_name not in window_snpid_to_reference_index:
			print('skip snp ' + gene_snpid_name)
			continue
		window_position, sign = window_snpid_to_reference_index[gene_snpid_name]

		if sign_aware_model:
			window_susie_mu[:, window_position] = (gene_susie_mu[:, ii])*sign
		else:
			window_susie_mu[:, window_position] = (gene_susie_mu[:, ii])

	return window_susie_mu

def get_window_regression_snp_indices(G_window_snpids, regression_snp_id_to_regression_snp_position):
	window_regression_snp_indices = []
	global_regression_snp_positions = []

	for snpid in G_window_snpids:
		if snpid in regression_snp_id_to_regression_snp_position:
			window_regression_snp_indices.append(True)
			global_regression_snp_positions.append(regression_snp_id_to_regression_snp_position[snpid])
		else:
			window_regression_snp_indices.append(False)

	# Put into nice array format
	window_regression_snp_indices = np.asarray(window_regression_snp_indices)
	global_regression_snp_positions = np.asarray(global_regression_snp_positions)

	# Quick error checking
	if np.sum(window_regression_snp_indices) != len(global_regression_snp_positions):
		print('assumption eroror')
		pdb.set_trace()

	return window_regression_snp_indices, global_regression_snp_positions

def compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld):
	gene_var = 0.0

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = (susie_mu)*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)


	num_susie_components = susie_mu.shape[0]
	for k_index in range(num_susie_components):
		gene_var = gene_var + np.sum((np.square(susie_mu[k_index,:]) + np.square(susie_mu_sd[k_index,:]))*np.diag(ld)*susie_alpha[k_index,:])
		eqtl_component_pmces = (susie_mu[k_index,:])*(susie_alpha[k_index,:])
		gene_var = gene_var - np.dot(np.dot(eqtl_component_pmces,ld), eqtl_component_pmces)
	gene_var = gene_var + np.dot(np.dot(gene_eqtl_effect_sizes,ld), gene_eqtl_effect_sizes)
				
	return gene_var


def extract_ld_annotations_for_this_gene_region(ld, susie_mu, susie_mu_sd, susie_alpha, variant_indices, n_ref_panel_samples):
	# Number of variants
	num_var = ld.shape[0]

	# Component level eqtl effect sizes for this gene		
	gene_component_effect_sizes = susie_mu*susie_alpha

	# eQTL effect sizes for this gene
	gene_eqtl_effect_sizes = np.sum(gene_component_effect_sizes,axis=0)

	# Compute squared eqtl effect sizes for this gene
	gene_squared_eqtl_effect_sizes = np.sum((np.square(susie_mu) + np.square(susie_mu_sd))*susie_alpha,axis=0) + gene_eqtl_effect_sizes*gene_eqtl_effect_sizes - np.sum(gene_component_effect_sizes*gene_component_effect_sizes,axis=0)

	# E[beta_k*beta_j]
	cross_terms = np.dot(np.reshape(gene_eqtl_effect_sizes, (num_var,1)), np.reshape(gene_eqtl_effect_sizes, (1,num_var))) - np.dot(np.transpose(gene_component_effect_sizes), gene_component_effect_sizes)

	# Compute gene variance
	gene_variance = compute_gene_variance(susie_mu, susie_mu_sd, susie_alpha, ld)

	# Compute ld scores for diagonal piece
	diagonal_ld_scores = np.sum((np.square(ld[variant_indices,:])*gene_squared_eqtl_effect_sizes),axis=1)
	
	# Comptute ld scores for off diagonal elements
	np.fill_diagonal(cross_terms, np.zeros(num_var))
	non_diagonal_ld_scores = np.sum(np.dot(cross_terms, ld[:, variant_indices])*ld[:,variant_indices],axis=0)

	# Generate complete ld scores
	ld_scores = (diagonal_ld_scores + non_diagonal_ld_scores)/gene_variance

	# Get adjusted ld scores
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))

	return ld_scores, adj_ld_scores

def extract_ld_annotations_for_this_gene_region_with_eqtl_point_estimate(ld, gene_eqtl_effect_sizes, variant_indices, n_ref_panel_samples, component_annotation):
	gene_variance = np.dot(np.dot(gene_eqtl_effect_sizes, ld), gene_eqtl_effect_sizes)

	standardized_gene_eqtl_effect_sizes = gene_eqtl_effect_sizes/np.sqrt(gene_variance)
	ld_scores = np.square(np.dot(ld[variant_indices,:], standardized_gene_eqtl_effect_sizes))

	# Get adjusted ld scores
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))

	anno_ld_scores = np.dot(np.reshape(ld_scores, (len(ld_scores),1)), np.reshape(component_annotation, (1, len(component_annotation))))
	anno_adj_ld_scores = np.dot(np.reshape(adj_ld_scores, (len(adj_ld_scores),1)), np.reshape(component_annotation, (1, len(component_annotation))))

	return anno_ld_scores, anno_adj_ld_scores

def extract_ld_annotations_for_this_gene_region_with_sample_correlation(pmces, sample_geno, variant_indices, n_ref_panel_samples):
	stand_sample_geno = np.copy(sample_geno)
	for kk in range(stand_sample_geno.shape[1]):
		stand_sample_geno[:, kk] = (sample_geno[:, kk] - np.mean(sample_geno[:,kk]))/np.std(sample_geno[:,kk])

	pred_expr = np.dot(stand_sample_geno, pmces)

	full = np.hstack((np.reshape(pred_expr, (len(pred_expr), 1)), stand_sample_geno[:, variant_indices]))

	tmper = np.corrcoef(np.transpose(full))

	ld_scores = np.square(tmper[0,1:])
	adj_ld_scores = ld_scores - ((1.0-ld_scores)/(n_ref_panel_samples-2.0))
	return ld_scores, adj_ld_scores

def extract_expected_component_annotation(snp_weights, G_window_snpids, snp_id_to_annotation_vec):
	# Error checking
	if len(snp_weights) != len(G_window_snpids):
		print('assumption erororo')
		pdb.set_trace()
	if np.abs(np.sum(snp_weights) - 1.0) > 1e-8:
		print('Less than one snp weights')
		print(np.sum(snp_weights))

	num_snps = len(snp_weights)
	for ii in range(num_snps):
		snp_id = G_window_snpids[ii]
		snp_anno = snp_id_to_annotation_vec[snp_id]
		if ii == 0:
			component_annotation = np.zeros(len(snp_anno))
		if snp_weights[ii] == 0.0:
			continue
		component_annotation = component_annotation + snp_anno*snp_weights[ii]
		
	return component_annotation

def extract_expected_common_component_annotation(snp_weights, G_window_snpids, snp_id_to_annotation_vec, snp_id_to_common_indicator):
	num_snps = len(snp_weights)
	for ii in range(num_snps):
		snp_id = G_window_snpids[ii]
		snp_anno = snp_id_to_annotation_vec[snp_id]
		common_var_indicator = snp_id_to_common_indicator[snp_id]
		if ii == 0:
			component_annotation = np.zeros(len(snp_anno))
		if snp_weights[ii] == 0.0:
			continue
		if common_var_indicator == 1.0:
			component_annotation = component_annotation + snp_anno*snp_weights[ii]
	return component_annotation



def compute_gene_level_ld_scores_for_single_gene(gene_tissue_model, ref_genotype_obj, snpid_to_reference_index, regression_snp_id_to_regression_snp_position, global_gene_level_ld_scores, global_gene_level_adj_ld_scores, snp_id_to_annotation_vec, snp_id_to_common_indicator, annotation_counter, annotation_counter_5_50):
	# First get gene model ub and lb snp name
	gene_snpid_names = np.asarray(gene_tissue_model['variant_names'])[:,0]
	gene_snp_positions = get_gene_model_snp_positions(gene_snpid_names)
	gene_snpid_lb = gene_snpid_names[np.argmin(gene_snp_positions)]
	gene_snpid_ub = gene_snpid_names[np.argmax(gene_snp_positions)]

	# Now get ref positions of gene snp id lb and ub
	ref_pos_gene_snp_lb = snpid_to_reference_index[gene_snpid_lb][0]
	ref_pos_gene_snp_ub = snpid_to_reference_index[gene_snpid_ub][0]

	# Gene window cm lb and ub
	gene_window_cm_lb = ref_genotype_obj['cm'][ref_pos_gene_snp_lb] - 1.0
	gene_window_cm_ub = ref_genotype_obj['cm'][ref_pos_gene_snp_ub] + 1.0

	# Get indices corresponding to ref genotype for this gene window
	ref_genotype_indices_for_window = (ref_genotype_obj['cm'] >= gene_window_cm_lb) & (ref_genotype_obj['cm'] < gene_window_cm_ub)

	# Subset ref genotype data to only snps in ref_genotype_indices_for_window
	G_window_geno = ref_genotype_obj['G'][:, ref_genotype_indices_for_window]
	G_window_ld = np.corrcoef(np.transpose(G_window_geno))
	G_window_rsids = ref_genotype_obj['rsid'][ref_genotype_indices_for_window]
	G_window_snpids = ref_genotype_obj['snp_id'][ref_genotype_indices_for_window]
	G_window_alt_snpids = ref_genotype_obj['alt_snp_id'][ref_genotype_indices_for_window]

	# Get number of window snps
	n_window_snps = len(G_window_alt_snpids)
	# Number of reference panel samples
	n_ref_panel_samples = G_window_geno.shape[0]

	# Get indices in window corresponding to regression snps
	window_regression_snp_indices, global_regression_snp_positions = get_window_regression_snp_indices(G_window_snpids, regression_snp_id_to_regression_snp_position)

	# Calculate number of regression snps in this window
	n_regression_snps_in_window = np.sum(window_regression_snp_indices)

	# If no regression snps, can return here
	if n_regression_snps_in_window == 0:
		print('SKIP WINDOW')
		return global_gene_level_ld_scores, global_gene_level_adj_ld_scores, annotation_counter, annotation_counter_5_50

	# Create mapping from snp id to reference snp position in this window
	window_snpid_to_reference_index = create_mapping_snpid_to_reference_index(G_window_snpids, G_window_alt_snpids)

	# Project gene models onto full window
	window_susie_mu = project_gene_model_onto_full_window(np.asarray(gene_tissue_model['susie_mu']), gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=True)
	window_susie_alpha = project_gene_model_onto_full_window(np.asarray(gene_tissue_model['susie_alpha']), gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=False)
	local_susie_mu_sd = np.sqrt(np.asarray(gene_tissue_model['susie_mu2']) - np.square(np.asarray(gene_tissue_model['susie_mu'])))
	window_susie_mu_sd = project_gene_model_onto_full_window(local_susie_mu_sd, gene_snpid_names, window_snpid_to_reference_index, n_window_snps, sign_aware_model=False)

	component_annotations = []
	component_annotations_5_50 = []

	component_bool = np.asarray(gene_tissue_model['component_bool'])[0,0]
	if component_bool:
		gene_components = np.asarray(gene_tissue_model['susie_cs'])[:,0] - 1
		for gene_component in gene_components:
			# Extract expected annotation vector for this meta snp
			snp_weights = window_susie_alpha[gene_component,:]
			component_annotation = extract_expected_component_annotation(snp_weights, G_window_snpids, snp_id_to_annotation_vec)
			component_annotation_5_50 = extract_expected_common_component_annotation(snp_weights, G_window_snpids, snp_id_to_annotation_vec, snp_id_to_common_indicator)

			component_annotations.append(component_annotation)
			component_annotations_5_50.append(component_annotation_5_50)

			# Extract LD scores using point estimate eqtl effect sizes
			component_pmces = window_susie_alpha[gene_component,:]*window_susie_mu[gene_component,:]
			window_gene_ld_scores, window_adj_ld_scores = extract_ld_annotations_for_this_gene_region_with_eqtl_point_estimate(G_window_ld, component_pmces, window_regression_snp_indices, n_ref_panel_samples, component_annotation)

			# Add updates to global array
			global_gene_level_ld_scores[global_regression_snp_positions, :] = global_gene_level_ld_scores[global_regression_snp_positions, :] + window_gene_ld_scores
			global_gene_level_adj_ld_scores[global_regression_snp_positions, :] = global_gene_level_adj_ld_scores[global_regression_snp_positions, :] + window_adj_ld_scores

			# Keep track of gene annotations
			annotation_counter = annotation_counter + component_annotation
			annotation_counter_5_50 = annotation_counter_5_50 + component_annotation_5_50

	return global_gene_level_ld_scores, global_gene_level_adj_ld_scores, annotation_counter, annotation_counter_5_50, component_annotations, component_annotations_5_50


def load_in_regression_snp_ids(variant_level_ld_score_file, rsid_to_snpid, snpid_to_reference_index):
	rsids = []
	snpids = []
	snp_to_regression_snp_position = {}
	f = gzip.open(variant_level_ld_score_file)
	head_count = 0
	snp_counter = 0
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		# This filter occurs because of weird filter at the beginnong on removing repeat snps (only allow this to happen once per genome)
		snp_id = rsid_to_snpid[rsid]
		snp_id_info = snp_id.split('_')
		snp_id_alt = snp_id_info[0] + '_' + snp_id_info[1] + '_' + snp_id_info[3] + '_' + snp_id_info[2] + '_' + snp_id_info[4]
		# Quick error check
		if snp_id not in snpid_to_reference_index:
			print('assumption eroror')
			pdb.set_trace()
	
		# Add to arrays
		rsids.append(rsid)
		snpids.append(snp_id)
		if snp_id in snp_to_regression_snp_position or snp_id_alt in snp_to_regression_snp_position:
			print('assumption eroror')
			pdb.set_trace()
		snp_to_regression_snp_position[snp_id] = snp_counter
		snp_to_regression_snp_position[snp_id_alt] = snp_counter
		snp_counter = snp_counter + 1
	f.close()

	return np.asarray(rsids), np.asarray(snpids), snp_to_regression_snp_position


def get_non_repeat_columns_from_plink_obj(G_obj, variant_level_ld_score_file):
	f = gzip.open(variant_level_ld_score_file)
	head_count = 0
	hm3_rs_ids = {}
	for line in f:
		line = line.decode('utf-8').rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rsid = data[1]
		if rsid in hm3_rs_ids:
			print('assumption erorro')
			pdb.set_trace()
		hm3_rs_ids[rsid] = 1
	f.close()

	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	G_obj_rs = np.asarray(G_obj.snp)

	n_var = len(G_obj_pos)
	used = {}
	valid_columns = []

	# Go through hm3 snps first
	for ii in range(n_var):
		snp_id1 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii]
		snp_id2 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a1[ii] + '_' + G_obj_a0[ii]
		rsid = G_obj_rs[ii]
		if rsid not in hm3_rs_ids:
			continue
		if snp_id1 in used or snp_id2 in used:
			print('assumption error')
			pdb.set_trace()
		used[snp_id1] = 1
		used[snp_id2] = 1
		valid_columns.append(ii)
	for ii in range(n_var):
		snp_id1 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a0[ii] + '_' + G_obj_a1[ii]
		snp_id2 = 'chr' + G_obj_chrom[ii] + '_' + str(G_obj_pos[ii]) + '_' + G_obj_a1[ii] + '_' + G_obj_a0[ii]
		rsid = G_obj_rs[ii]
		if rsid in hm3_rs_ids:
			continue
		if snp_id1 in used or snp_id2 in used:
			continue
		used[snp_id1] = 1
		used[snp_id2] = 1
		valid_columns.append(ii)
	valid_columns = np.sort(np.asarray(valid_columns))
	return valid_columns

def create_mapping_from_snpid_to_annotation_vec(variant_level_annotation_file, rsid_to_snpid):
	f = open(variant_level_annotation_file)
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			anno_names = np.asarray(data[4:])
			continue
		rsid = data[2]
		if rsid not in rsid_to_snpid:
			continue
		annotation_vec = np.asarray(data[4:]).astype(float)
		snp_id = rsid_to_snpid[rsid]
		snp_info = snp_id.split('_')
		snp_id_alt = snp_info[0] + '_' + snp_info[1] + '_' + snp_info[3] + '_' + snp_info[2] + '_' + snp_info[4]

		if snp_id in dicti or snp_id_alt in dicti:
			print('assumption eroror!')
			pdb.set_trace()
		dicti[snp_id] = annotation_vec
		dicti[snp_id_alt] = annotation_vec
	return dicti, anno_names

def compute_mafs_of_reference_snps(G_obj_geno):
	mafs = []
	n_snps = G_obj_geno.shape[1]
	n_ind = G_obj_geno.shape[0]
	for snp_iter in range(n_snps):
		af = np.sum(G_obj_geno[:, snp_iter])/(2.0*n_ind)
		if af > .5:
			maf = 1.0 - af 
		else:
			maf = af
		mafs.append(maf)
	return np.asarray(mafs)

def create_mapping_from_snp_id_to_common_boolean(G_obj_snp_ids, G_obj_alt_snp_ids, G_obj_mafs):
	n_snps = len(G_obj_snp_ids)
	dicti = {}
	for ii in range(n_snps):
		snp_id = G_obj_snp_ids[ii]
		alt_snp_id = G_obj_alt_snp_ids[ii]
		maf = G_obj_mafs[ii]
		booler = 0.0
		if maf > .05:
			booler = 1.0
		dicti[snp_id] = booler
		dicti[alt_snp_id] = booler
	return dicti


variant_level_annotation_file = sys.argv[1]
variant_level_ld_score_file = sys.argv[2]
genotype_stem = sys.argv[3]
pseudotissue_name = sys.argv[4]
gtex_susie_gene_models_dir = sys.argv[5]
chrom_num = sys.argv[6]
gene_type = sys.argv[7]
gene_level_sldsc_output_root = sys.argv[8]


print('CHROM' + str(chrom_num))


# Create dictionary from pseudotissue to array of gene models for that tissue
tissue_to_gene_model_df = create_tissue_to_gene_model_df(pseudotissue_name, gtex_susie_gene_models_dir, gene_type, chrom_num)

# Number of genes on this chromosome in this tissue
n_genes = tissue_to_gene_model_df.shape[0]

# Load in Reference Genotype data
G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
# For some reason, this file has some repeats. First get columns that don't correspond to repeats
valid_variant_columns = get_non_repeat_columns_from_plink_obj(G_obj, variant_level_ld_score_file)
G_obj = G_obj[:, valid_variant_columns]

G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
G_obj_chrom = np.asarray(G_obj.chrom)
G_obj_pos = np.asarray(G_obj.pos)
# For our purposes, a0 is the effect allele
# For case of plink package, a0 is the first column in the plink bim file
G_obj_a0 = np.asarray(G_obj.a0)
G_obj_a1 = np.asarray(G_obj.a1)
# RSids
G_obj_rsids = np.asarray(G_obj.snp)
# Centimorgan distances
G_obj_cm = np.asarray(G_obj.cm)
# Snp ids
G_obj_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a0 + '_' + G_obj_a1 + '_b38'
G_obj_alt_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a1 + '_' + G_obj_a0 + '_b38'

# Compute mafs of of reference snps
G_obj_mafs = compute_mafs_of_reference_snps(G_obj_geno)

# Create mapping from snpid to common
snp_id_to_common_indicator = create_mapping_from_snp_id_to_common_boolean(G_obj_snp_ids, G_obj_alt_snp_ids, G_obj_mafs)

# Put geno into organized dictionary
genotype_obj = {'G': G_obj_geno, 'rsid': G_obj_rsids, 'snp_id': G_obj_snp_ids, 'alt_snp_id': G_obj_alt_snp_ids, 'position': G_obj_pos, 'cm': G_obj_cm}

# Quick error checks on genotype data
if len(np.unique(G_obj_rsids)) != len(G_obj_rsids):
	print('assumption eroror')
	pdb.set_trace()
if len(np.unique(G_obj_snp_ids)) != len(G_obj_snp_ids):
	print('assumption error')
	pdb.set_trace()


# Create mapping from rsid to snpid and snpid to rsid
rsid_to_snpid = create_mapping_from_rsid_to_snpid(genotype_obj['rsid'], genotype_obj['snp_id'])
snpid_to_rsid = create_mapping_from_snpid_to_rsid(genotype_obj['rsid'], genotype_obj['snp_id'])

# Create mapping from SNPID to annotation vec
snp_id_to_annotation_vec, anno_names = create_mapping_from_snpid_to_annotation_vec(variant_level_annotation_file, rsid_to_snpid)

# Create mapping from snp id to reference index
snpid_to_reference_index = create_mapping_snpid_to_reference_index(genotype_obj['snp_id'], genotype_obj['alt_snp_id'])

# Get regression rsids
regression_rsids, regression_snpids, regression_snp_id_to_regression_snp_position = load_in_regression_snp_ids(variant_level_ld_score_file, rsid_to_snpid, snpid_to_reference_index)
n_regression_snps = len(regression_rsids)

# Initialize vector to keep track of gene level ld scores for regression snps on this chromosome
global_gene_level_ld_scores = np.zeros((n_regression_snps, len(anno_names)))
global_gene_level_adj_ld_scores = np.zeros((n_regression_snps, len(anno_names)))

# Initialize vector to keep track of annotations
annotation_counter = np.zeros(len(anno_names))
annotation_counter_5_50 = np.zeros(len(anno_names))

# Initialize file to keep track of annotations
annotation_output_file = gene_level_sldsc_output_root + pseudotissue_name + '_gene_annotations.txt'
t = open(annotation_output_file,'w')
annotation_5_50_output_file = gene_level_sldsc_output_root + pseudotissue_name + '_gene_annotations_5_50.txt'
t2 = open(annotation_5_50_output_file,'w')

start_time = time.time()
# Loop through genes
for g_counter in range(n_genes):
	#print(g_counter)
	gene_info_vec = tissue_to_gene_model_df[g_counter,:]
	rdat_gene_model_file = gene_info_vec[0]
	# Load gene-tissue model from gene-tissue weight file
	gene_tissue_model = pyreadr.read_r(rdat_gene_model_file)

	global_gene_level_ld_scores, global_gene_level_adj_ld_scores, annotation_counter, annotation_counter_5_50, component_annotations, component_annotations_5_50 = compute_gene_level_ld_scores_for_single_gene(gene_tissue_model, genotype_obj, snpid_to_reference_index, regression_snp_id_to_regression_snp_position, global_gene_level_ld_scores, global_gene_level_adj_ld_scores, snp_id_to_annotation_vec, snp_id_to_common_indicator, annotation_counter, annotation_counter_5_50)

	end_time = time.time()
	#print(end_time-start_time)
	gene_id = gene_info_vec[1]

	for comp_iter, gene_anno_vec in enumerate(component_annotations):
		t.write(gene_id + '_' + str(comp_iter) + '\t' + '\t'.join(gene_anno_vec.astype(str)) + '\n')
	for comp_iter, gene_anno_vec in enumerate(component_annotations_5_50):
		t2.write(gene_id + '_' + str(comp_iter) + '\t' + '\t'.join(gene_anno_vec.astype(str)) + '\n')

	start_time = end_time

t.close()
t2.close()

# Save results to output file
output_file1 = gene_level_sldsc_output_root + pseudotissue_name + '_gene_ld_scores'
np.savetxt(output_file1, global_gene_level_ld_scores, fmt="%s", delimiter='\t')

output_file2 = gene_level_sldsc_output_root + pseudotissue_name + '_gene_adj_ld_scores'
np.savetxt(output_file2, global_gene_level_adj_ld_scores, fmt="%s", delimiter='\t')


output_file3 = gene_level_sldsc_output_root + pseudotissue_name + '_gene_m'
np.savetxt(output_file3, annotation_counter, fmt="%s", delimiter='\t')

output_file4 = gene_level_sldsc_output_root + pseudotissue_name + '_gene_m_5_50'
np.savetxt(output_file4, annotation_counter_5_50, fmt="%s", delimiter='\t')




