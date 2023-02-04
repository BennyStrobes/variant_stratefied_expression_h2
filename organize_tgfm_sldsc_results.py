import numpy as np 
import os
import sys
import pdb



def get_anno_names(example_anno_file):
	f = open(example_anno_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		anno_names = np.asarray(data[3:])
		break
	f.close()
	return anno_names




def extract_tau_and_tau_se_from_log_file(sldsc_log_file, anno_names):
	f = open(sldsc_log_file)
	for line in f:
		line = line.rstrip()
		if line.startswith('Coefficients:'):
			coef = np.asarray(line.split()[1:]).astype(float)
			if len(coef) < len(anno_names):
				line = f.next()
				data = np.asarray(line.split()).astype(float)
				coef = np.hstack([coef, data])
				if len(coef) < len(anno_names):
					line = f.next()
					data = np.asarray(line.split()).astype(float)
					coef = np.hstack([coef, data])				
		if line.startswith('Coefficient SE:'):
			coef_se = np.asarray(line.split()[2:]).astype(float)
			if len(coef_se) < len(anno_names):
				line = f.next()
				data = np.asarray(line.split()).astype(float)
				coef_se = np.hstack([coef_se, data])
				if len(coef_se) < len(anno_names):
					line = f.next()
					data = np.asarray(line.split()).astype(float)
					coef_se = np.hstack([coef_se, data])				
	return coef, coef_se

def extract_tau_and_tau_se_from_summary_file(sldsc_summary_file):
	f = open(sldsc_summary_file)
	head_count = 0
	tau = []
	tau_se = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tau.append(float(data[-3]))
		tau_se.append(float(data[-2]))
	f.close()
	return np.asarray(tau), np.asarray(tau_se)


def get_anno_counts(anno_stem):
	for chrom_num in range(1,23):
		filer = anno_stem + str(chrom_num) + '.l2.M_5_50'
		counts = np.loadtxt(filer)
		if chrom_num == 1:
			totaler = counts
		else:
			totaler = totaler + counts

	return totaler

def get_binary_annotations(anno_stem):
	binary = np.asarray([True]*93)
	for chrom_num in range(1,3):
		anno_file = anno_stem + str(chrom_num) + '.annot'
		f = open(anno_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			anno = data[4:]
			variant_anno = np.asarray(anno[:93]).astype(float)
			indices =(variant_anno!=1) & (variant_anno!=0)
			binary = (~indices)*(binary)
		f.close()
	return binary

def extract_annotation_sdevs(anno_stem):
	var_anno_arr = []
	gene_anno_arr = []
	for chrom_num in range(1,23):
		anno_file = anno_stem + str(chrom_num) + '.annot'
		f = open(anno_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			anno = data[4:]
			if data[2].startswith('ENSG') == True:
				gene_anno = np.asarray(anno[94:]).astype(float)
				gene_anno_arr.append(gene_anno)
			else:
				variant_anno = np.asarray(anno[1:93]).astype(float)
				var_anno_arr.append(variant_anno)
		f.close()
	var_anno_arr = np.asarray(var_anno_arr)
	gene_anno_arr = np.asarray(gene_anno_arr)
	var_sdev = np.std(var_anno_arr,axis=0)
	gene_sdev = np.std(gene_anno_arr,axis=0)
	return var_sdev, gene_sdev



sldsc_output_root = sys.argv[1]
anno_stem = sys.argv[2]


sldsc_log_file = sldsc_output_root + '.log'
example_anno_file = anno_stem + '21.l2.ldscore'
sldsc_summary_file = sldsc_output_root + '.results'

# Extract annotation names
anno_names = get_anno_names(example_anno_file)

# Extract taus and z-scores
tau, tau_se = extract_tau_and_tau_se_from_summary_file(sldsc_summary_file)
tau_z = tau/tau_se

# Compute mediated h2 for each annotation
anno_counts = get_anno_counts(anno_stem)
t = open(sldsc_output_root + 'organized_mediated_h2.txt','w')
t.write('Annotation\tmediated_h2\tmediated_h2_se\tmediated_h2_z\n')
partitioned_h2 = tau*anno_counts
partitioned_h2_se = tau_se*anno_counts
partitioned_h2_z = partitioned_h2/partitioned_h2_se
for ii in range(len(anno_names)):
	t.write(anno_names[ii] + '\t' + str(partitioned_h2[ii]) + '\t' + str(partitioned_h2_se[ii]) + '\t' + str(partitioned_h2_z[ii]) + '\n')
t.close()
print(np.sum(partitioned_h2[93:])/np.sum(partitioned_h2))




# Compute tau* for all annotations
# Exctract sdevs of each annotation
variant_anno_sdevs, gene_anno_sdevs = extract_annotation_sdevs(anno_stem)
variant_anno_names = anno_names[1:93]
gene_anno_names = anno_names[94:]
variant_taus = tau[1:93]
gene_taus = tau[94:]
variant_tau_ses = tau_se[1:93]
gene_tau_ses = tau_se[94:]
n_var_anno = len(variant_anno_sdevs)
n_gene_anno = len(gene_anno_names)
t = open(sldsc_output_root + 'tau_star.txt','w')
t.write('annotation_name\ttau_star\ttau_star_se\ttau_star_z\n')
for var_anno_iter in range(n_var_anno):
	per_variant_h2 = np.sum(partitioned_h2[:93])/anno_counts[0]
	scaling_factor = variant_anno_sdevs[var_anno_iter]/per_variant_h2
	tau_star = variant_taus[var_anno_iter]*scaling_factor
	tau_star_se = variant_tau_ses[var_anno_iter]*scaling_factor
	tau_star_z = tau_star/tau_star_se
	t.write(variant_anno_names[var_anno_iter] + '\t' + str(tau_star) + '\t' + str(tau_star_se) + '\t' + str(tau_star_z) + '\n')
for gene_anno_iter in range(n_gene_anno):
	per_gene_h2 = np.sum(partitioned_h2[93:])/anno_counts[93]
	scaling_factor = gene_anno_sdevs[gene_anno_iter]/per_gene_h2
	tau_star = gene_taus[gene_anno_iter]*scaling_factor
	tau_star_se = gene_tau_ses[gene_anno_iter]*scaling_factor
	tau_star_z = tau_star/tau_star_se
	t.write(gene_anno_names[gene_anno_iter] + '\t' + str(tau_star) + '\t' + str(tau_star_se) + '\t' + str(tau_star_z) + '\n')
t.close()



# Print enrichments for exclusively binary annotations
binary_annotations = get_binary_annotations(anno_stem)
t = open(sldsc_output_root + 'organized_binary_enrichments.txt','w')
t.write('Annotation\tprop_snps\tenrichment\tenrichment_se\tenrichment_p\n')
f = open(sldsc_summary_file)
head_count = 0
line_counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	mod_line_counter = np.mod(line_counter, 93)
	line_counter = line_counter + 1
	if binary_annotations[mod_line_counter] == False:
		continue
	if mod_line_counter == 0:
		continue
	anno_name = data[0]
	prop_snps = data[1]
	enrichment = data[4]
	enrichment_se = data[5]
	enrichment_p = data[6]
	t.write(anno_name + '\t' + prop_snps + '\t' + enrichment + '\t' + enrichment_se + '\t' + enrichment_p + '\n')
f.close()
t.close()


