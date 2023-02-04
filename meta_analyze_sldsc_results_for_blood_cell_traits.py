import numpy as np 
import os
import sys
import pdb


def load_in_tau_star_data_across_traits(trait_names, sldsc_results_dir):
	tau_star_arr = []
	tau_star_se_arr = []
	anno_names_arr = []
	for trait_name in trait_names:
		tau_star_file = sldsc_results_dir + trait_name + '_Whole_Blood_baselineLD_no_qtl_gene_adj_ld_scores_tau_star.txt'
		data = np.loadtxt(tau_star_file,dtype=str,delimiter='\t')
		anno_names = data[1:,0]
		tau_star = data[1:,1].astype(float)
		tau_star_se = data[1:,2].astype(float)
		tau_star_arr.append(tau_star)
		tau_star_se_arr.append(tau_star_se)
		anno_names_arr.append(anno_names)
	return np.transpose(np.asarray(tau_star_arr)), np.transpose(np.asarray(tau_star_se_arr)), anno_names_arr[0]

def load_in_enrichment_data_across_traits(trait_names, sldsc_results_dir):
	tau_star_arr = []
	tau_star_se_arr = []
	anno_names_arr = []
	for trait_name in trait_names:
		tau_star_file = sldsc_results_dir + trait_name + '_Whole_Blood_baselineLD_no_qtl_gene_adj_ld_scores_organized_binary_enrichments.txt'
		data = np.loadtxt(tau_star_file,dtype=str,delimiter='\t')
		anno_names = data[1:,0]
		tau_star = data[1:,2].astype(float)
		tau_star_se = data[1:,3].astype(float)
		tau_star_arr.append(tau_star)
		tau_star_se_arr.append(tau_star_se)
		anno_names_arr.append(anno_names)
	return np.transpose(np.asarray(tau_star_arr)), np.transpose(np.asarray(tau_star_se_arr)), anno_names_arr[0]


def load_in_sldsc_results_data_across_traits(trait_names, sldsc_results_dir):
	# Load in tau_star data across traits
	tau_star_mat, tau_star_se_mat, tau_star_anno_names = load_in_tau_star_data_across_traits(trait_names, sldsc_results_dir)
	# Load in enrichment data across traits
	enrichment_mat, enrichment_se_mat, enrichment_anno_names = load_in_enrichment_data_across_traits(trait_names, sldsc_results_dir)

	return tau_star_mat, tau_star_se_mat, enrichment_mat, enrichment_se_mat, tau_star_anno_names, enrichment_anno_names

def meta_analysis(effects, se, method='random', weights=None):
	# From Omer Weissbrod
	assert method in ['fixed', 'random']
	d = effects
	variances = se**2

	#compute random-effects variance tau2
	vwts = 1.0 / variances
	fixedsumm = vwts.dot(d) / vwts.sum()
	Q = np.sum(((d - fixedsumm)**2) / variances)
	df = len(d)-1
	tau2 = np.maximum(0, (Q-df) / (vwts.sum() - vwts.dot(vwts) / vwts.sum()))

	#defing weights
	if weights is None:
		if method == 'fixed':
			wt = 1.0 / variances
		else:
			wt = 1.0 / (variances + tau2)
	else:
		wt = weights

	#compute summtest
	summ = wt.dot(d) / wt.sum()
	if method == 'fixed':
		varsum = np.sum(wt*wt*variances) / (np.sum(wt)**2)
	else:
		varsum = np.sum(wt*wt*(variances+tau2)) / (np.sum(wt)**2)
	###summtest = summ / np.sqrt(varsum)

	summary=summ
	se_summary=np.sqrt(varsum)

	return summary, se_summary

def run_meta_analysis(tau_star_mat, tau_star_se_mat, tau_star_anno_names, meta_analyzed_tau_star_output, value_name):
	t = open(meta_analyzed_tau_star_output,'w')
	t.write('annotation_name\t' + value_name + '\t' + value_name + '_se\t' + value_name + '_z\n')

	for anno_iter, anno_name in enumerate(tau_star_anno_names):
		ma_tau_star, ma_tau_star_std_error = meta_analysis(tau_star_mat[anno_iter,:], tau_star_se_mat[anno_iter,:], method='random')
		ma_tau_star_z = ma_tau_star/ma_tau_star_std_error

		t.write(anno_name + '\t' + str(ma_tau_star) + '\t' + str(ma_tau_star_std_error) + '\t' + str(ma_tau_star_z) + '\n')

	t.close()



sldsc_results_dir = sys.argv[1]


trait_names = ['blood_RBC_DISTRIB_WIDTH', 'blood_RED_COUNT', 'blood_WHITE_COUNT', 'blood_PLATELET_COUNT', 'blood_EOSINOPHIL_COUNT']


# Load in cross trait data
tau_star_mat, tau_star_se_mat, enrichment_mat, enrichment_se_mat, tau_star_anno_names, enrichment_anno_names = load_in_sldsc_results_data_across_traits(trait_names, sldsc_results_dir)

# Meta analyze tau stars
meta_analyzed_tau_star_output = sldsc_results_dir + 'blood_trait_meta_analysis_Whole_Blood_baselineLD_no_qtl_gene_adj_ld_scores_tau_star.txt'
run_meta_analysis(tau_star_mat, tau_star_se_mat, tau_star_anno_names, meta_analyzed_tau_star_output, 'tau_star')


# Meta analyze enrichments
meta_analyzed_enrichment_output = sldsc_results_dir + 'blood_trait_meta_analysis_Whole_Blood_baselineLD_no_qtl_gene_adj_ld_scores_binary_enrichments.txt'
run_meta_analysis(enrichment_mat, enrichment_se_mat, enrichment_anno_names, meta_analyzed_enrichment_output, 'enrichment')


