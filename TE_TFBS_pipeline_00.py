import os
import subprocess
import pandas as pd
from pandas import DataFrame
import itertools
import collections
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from numpy import sqrt, abs, round
from scipy.stats import norm
import sys
import argparse

###############################################
###This script is used to investigate the regulatory activity of a transposon subfamily in specific diseases or developmental stage.
###This script requires input files including a transcription factor binding sites (TFBS) bed file (bed4), a TE bed file (bed6), and the path to the directory where the outputs will be saved to
###This script will do the following:
###1)Identify TFBSs that overlap with TEs of interest
###2)Identity the number of occurrences of TFBS colocalisation
###3)Calculate the enrichment of TFBS colocalisation in TEs compared to the remainder of the genome
###4)Identify transcription start sites (TSS) in TF-bound TEs.

###############################################


###############################################
###Download the TFBS file data of interest, and investigate the overlaps between these TFBSs and the transposons of interest. 

# TFBS_bed_file = '/data/clara/L1PA2_pipeline_script/test_files/test_all_TFBS_annotated.bed'
# TE_bed_file = '/data/clara/L1PA2_pipeline_script/test_files/test_L1PA2.bed'
# TEname = 'test_L1PA2'
# Output_dir = '/data/clara/L1PA2_pipeline_script/'
# TSS_bed_file = '/data/clara/PhD/L1PA2-derived_oncogenes/L1PA2_MCF7_MCF10A_RNA-seq/GSE108541_MCF7_vs_MCF10A_transcript_FPKM_Upreg_TSS.bed'
# TSS_name = 'Upreg'

def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Investigating the regulatory activity of TEs')
    parser.add_argument('-a', '--TF_bed', action="store", help='absolute path to TFBS bed file (bed4)')
    parser.add_argument('-b', '--TE_bed', action="store", help='absolute path to TE bed file (bed6)')
    parser.add_argument('-t', '--TE_name', action="store", help='name of the TE subfmaily')
    parser.add_argument('-o', '--out_dir', action="store", help='absolute path for the output directory')
    parser.add_argument('-r', '--TSS_bed', action="store", help='absolute path for the TSS bed file')
    parser.add_argument('-l', '--TSS_name', action="store", help='label for the TSS bed file')
    options = parser.parse_args(args)
    return options

options = getOptions(sys.argv[1:])

###############################################
###Checking for input files
TFBS_bed_file = options.TF_bed
TE_bed_file = options.TE_bed
Output_dir = options.out_dir
TEname = options.TE_name
TSS_bed_file = options.TSS_bed
TSS_label = options.TSS_name

if not TFBS_bed_file:
	sys.exit('Missing TFBS file. Exiting...')
if not TE_bed_file:
	sys.exit('Missing TE file. Exiting...')
if not Output_dir:
	sys.exit('Missing output dir. Exiting...')
if not TEname:
	sys.exit('Missing TE name. Exiting...')
if not TSS_bed_file:
	sys.exit('Missing TSS bed file. Exiting...')
if not TSS_label:
	sys.exit('Missing TSS label. Exiting...')	

###############################################
###Checking bed file formats
def Counting_lines_function(file):
	out = subprocess.Popen(['wc', '-l', file], 
			stdout=subprocess.PIPE, 
			stderr=subprocess.STDOUT)
	stdout,stderr = out.communicate()
	count = stdout.split()[0].decode('ascii')	
	return count

if Counting_lines_function(TE_bed_file) == '0':
	sys.exit('The TE bed file is empty. Exiting...')
if Counting_lines_function(TFBS_bed_file) == '0':
	sys.exit('The TFBS bed file is empty. Exiting...')
if Counting_lines_function(TSS_bed_file) == '0':
	sys.exit('The TSS bed file is empty. Exiting...')	

def Checking_column_number(file):
	df = pd.read_csv(file,sep = '\t', header = None, index_col = None)
	column_number = len(df.columns)
	return column_number

if Checking_column_number(TE_bed_file) != 6:
	sys.exit('The TE bed file is not in a bed6 format. Exiting...')
if Checking_column_number(TFBS_bed_file) != 4:
	sys.exit('The TFBS bed file is not in a bed4 format. Exiting...')

###############################################
###Perform TE and TFBS intersection
Work_dir = Output_dir + TEname + '_TFBS_dir/'
dir_check = os.path.isdir(Work_dir)
if dir_check == True:
	sys.exit('The working directory already exists. Exiting...')

cmd = 'mkdir ' + Work_dir
os.system(cmd)

TFBS_int_dir = Work_dir + TEname + '_TFBS_intersection/'
cmd = 'mkdir ' + TFBS_int_dir
os.system(cmd)

print('Performing TFBS intersection with', TEname)

TE_TFBS_wa_file = TFBS_int_dir + TEname+'_TFBS_F0.5_wa.bed'
cmd = 'bedtools intersect -u -F 0.5 -wa -a '+TE_bed_file+ ' -b '+TFBS_bed_file + ' > ' + TE_TFBS_wa_file
os.system(cmd)

TE_TFBS_wb_file = TFBS_int_dir + TEname+'_TFBS_F0.5_wb.bed'
cmd = 'bedtools intersect -u -f 0.5 -wa -b '+TE_bed_file+ ' -a '+TFBS_bed_file + ' > ' + TE_TFBS_wb_file
os.system(cmd)	

TE_TFBS_wo_file = TFBS_int_dir + TEname+'_TFBS_F0.5_wo.txt'
cmd = 'bedtools intersect -F 0.5 -wo -a '+TE_bed_file+ ' -b '+TFBS_bed_file + ' > ' + TE_TFBS_wo_file
os.system(cmd)

print('Total_number of TEs in this subfamily', Counting_lines_function(TE_bed_file))
print('Number of TEs with at least one TFBS binding site', Counting_lines_function(TE_TFBS_wa_file))
if Counting_lines_function(TE_TFBS_wa_file) == '0':
	sys.exit('No TE overlaps with TFBSs. Exiting...')
print('Number of TFBSs that overlap with TEs', Counting_lines_function(TE_TFBS_wb_file))
if Counting_lines_function(TE_TFBS_wa_file) == '0':
	sys.exit('No TFBS overlaps with TEs. Exiting...')

###############################################
###Summarising the number of TEs that contained the binding sites of each TF
print('Summarising the number of TEs that contained the binding sites of each TF')
TE_TFBS_wo_df = pd.read_csv(TE_TFBS_wo_file, sep = '\t', names = ['TE_chromosome', 'TE_geno_start', 'TE_geno_end','TE_name', 'TE_Score', 'TE_strand', 'TF_chromosome', 'TF_geno_start', 'TF_geno_end','TF_name', 'Overlap'])

Unique_wo_df = TE_TFBS_wo_df.loc[:,['TE_chromosome', 'TE_geno_start', 'TE_geno_end','TF_name']]
Unique_wo_df['TE'] = Unique_wo_df['TE_chromosome'].map(str) + '-' + Unique_wo_df['TE_geno_start'].map(str) + ':' + Unique_wo_df['TE_geno_end'].map(str)
Unique_wo_df = Unique_wo_df.loc[:,['TE','TF_name']]
Unique_wo_df = Unique_wo_df.drop_duplicates()
TF_TE_count_df = Unique_wo_df['TF_name'].value_counts().to_frame().reset_index()
TF_TE_count_df.columns = ['TF_name', 'Number of bound TEs']

TF_TE_count_file = TFBS_int_dir + TEname+'_TFBS_occurrences_summary.bed'
TF_TE_count_df.to_csv(TF_TE_count_file, sep = '\t', header = True, index = False)


###############################################
###Investigate the pattern of TF binding co-localisation in these TEs
print('Investigating the pattern of TF binding co-localisation in these TEs')
TFBS_coloc_dir = Work_dir + TEname + '_TFBS_coloc/'
cmd = 'mkdir ' + TFBS_coloc_dir
os.system(cmd)

TE_TFBS_wb_df = pd.read_csv(TE_TFBS_wb_file, sep = '\t', names = ['chromosome', 'geno_start', 'geno_end','TF_name'])
TF_name_array = TE_TFBS_wb_df['TF_name'].to_numpy()
TF_name_list = list(set(TF_name_array))
TF_count = len(TF_name_list)
print('Number of TFs that bind to',TEname, TF_count)
print(TF_name_list)

###Count the occurrences of TFBS colocalisation
TF_pair_count_dict = {}
TF_pairs_list = list(itertools.product(TF_name_list,repeat=2))
for TF_pair in TF_pairs_list:
	TF_pair_count_dict[TF_pair] = 0

TE_TF_list_dict = dict()
with open (TE_TFBS_wo_file, 'r') as a:
	lines = a.readlines()
	for line in lines:
		line = line.strip()
		fields = line.split()
		
		TE_chromo = fields[0]
		TE_start = fields[1]
		TE_end = fields[2]
		TE_name = fields[3]
		TE_score = fields[4]
		TE_strand = fields[5]
		TE_key = '\t'.join([TE_chromo, TE_start, TE_end, TE_name, TE_score, TE_strand])
		TF_name = fields[9]

		if TE_key not in TE_TF_list_dict:
			TE_TF_list_dict[TE_key] = [TF_name]
		else:
			TE_TF_list_dict[TE_key].append(TF_name)

for each_key in TE_TF_list_dict.keys():
	TF_list = TE_TF_list_dict[each_key]
	TF_list = list(set(TF_list))
	TF_obs_pairs_list = list(itertools.product(TF_list,repeat=2))
	for TF_obs_pair in TF_obs_pairs_list:
		TF_pair_count_dict[TF_obs_pair] = TF_pair_count_dict[TF_obs_pair] + 1

count_list = []
permutation_outfile = TFBS_coloc_dir + TEname + '_TFBS_permutations_count.txt'
with open (permutation_outfile, 'w') as a:
	for key,value in TF_pair_count_dict.items():
		a.write(str(key)+'\t'+str(value)+'\n')
		count_list.append(value)

count_chunk_list = [count_list[x:x+TF_count] for x in range(0, len(count_list), TF_count)]
pair_df = DataFrame (count_chunk_list,  index = TF_name_list, columns = TF_name_list)

##Reorder the TF_names by the binding frequency
TF_binding_frequency_dict = {}
for TF_name in TF_name_list:
	key_list_temp = []
	key_list_temp.append(TF_name)
	key_pair_list_temp = list(itertools.product(key_list_temp,repeat=2))
	for key_pair in key_pair_list_temp:
		TF_binding_frequency_dict[TF_name] = TF_pair_count_dict[key_pair]

sorted_x = sorted(TF_binding_frequency_dict.items(), key=lambda kv: kv[1], reverse = True)
TF_binding_frequency_dict = collections.OrderedDict(sorted_x)

TF_key_ordered_list = []
for TF_name in TF_binding_frequency_dict:
	TF_key_ordered_list.append(TF_name)

df_ordered = pair_df.reindex(columns=TF_key_ordered_list)
df_ordered = df_ordered.reindex(index=TF_key_ordered_list)

dataframe_raw = TFBS_coloc_dir + TEname + '_TFBS_permutations_count_dataframe_raw.txt'
with open(dataframe_raw, 'a') as f:
	f.write(df_ordered[TF_key_ordered_list].to_string(header = TF_key_ordered_list, index = TF_key_ordered_list))

###Adjust each row by total number
df_adj_row = df_ordered
for TF_name in TF_key_ordered_list:
	df_adj_row.loc[TF_name] = df_adj_row.loc[TF_name]/TF_binding_frequency_dict[TF_name]

dataframe_adj_row = TFBS_coloc_dir + TEname + '_TFBS_permutations_count_dataframe_adj_row.txt'
with open(dataframe_adj_row, 'a') as f:
	f.write(df_adj_row[TF_key_ordered_list].to_string(header = TF_key_ordered_list, index = TF_key_ordered_list))

###Calculating the distances between the TFBSs
distance_list = []
with open (TE_TFBS_wo_file, 'r') as a:
	lines = a.readlines()
	for line in lines:
		line = line.strip()
		fields = line.split()

		TE_chromo = fields[0]
		TE_start = fields[1]
		TE_end = fields[2]
		TE_name = fields[3]
		TE_score = fields[4]
		TE_strand = fields[5]

		TF_start = int(fields[7])
		TF_end = int(fields[8])
		TF_name = fields[9]
		mid_TF_pos = int((TF_end + TF_start)/2)

		if TE_strand == '+':
			distance = mid_TF_pos - int(TE_start)
		else:
			distance = int(TE_end) - mid_TF_pos
		newline = '\t'.join([TE_chromo, TE_start, TE_end, TE_name, TE_score, TE_strand, TF_name, str(distance)])
		distance_list.append(newline)

distance_file =  TFBS_int_dir + TEname + '_TFBS_wo_F0.5_distance.txt'
with open (distance_file, 'w') as b:
	for line in distance_list:
		b.write(line + '\n')

TE_key_list = []
distance_dict = {}
with open (distance_file, 'r') as c:
	lines = c.readlines()
	for line in lines:
		line = line.strip()
		fields = line.split()

		TE_chromo = fields[0]
		TE_start = fields[1]
		TE_end = fields[2]
		TE_name = fields[3]
		TE_score = fields[4]
		TE_strand = fields[5]
		TE_key = '\t'.join([TE_chromo, TE_start, TE_end, TE_name, TE_score, TE_strand])
		TE_key_list.append(TE_key)

		TF_name = fields[6]
		TF_dist = fields[7]

		combined_key = TE_key+'|'+TF_name
		distance_dict[combined_key] = TF_dist

TE_key_list = list(set(TE_key_list))

def calculating_pair_distance(TF1, TF2):
	for TE in TE_key_list:
		try:
			combined_key_1 = TE+'|'+TF1
			combined_key_2 = TE+'|'+TF2

			pos1 = int(distance_dict[combined_key_1])
			pos2 = int(distance_dict[combined_key_2])
			pair_dist = abs(pos1 - pos2)

			Pair_distance_line = '\t'.join([TE, TF1,TF2, str(pair_dist)])
			Pair_distance_list.append(Pair_distance_line)

		except:
			KeyError
	return

Pair_distance_list = []
TF_pairs_list = list(itertools.combinations(TF_name_list, 2))
for TF_pair in TF_pairs_list:
	TF1 = TF_pair[0]
	TF2 = TF_pair[1]
	calculating_pair_distance(TF1, TF2)

Pair_distance_file = TFBS_coloc_dir + TEname + '_TFBS_permutations_distances.txt'
with open (Pair_distance_file, 'w') as a:
	for line in Pair_distance_list:
		a.write(line + '\n')

###############################################
###Evaluate the enrichment of TF binding co-localisation in these TEs
###Firstly, get the global genomic locations file for the TE-binding TFs
print('Evaluating the enrichment of TF binding co-localisation in these TEs')
TFBS_files_dir = Work_dir + TEname + '_TFBS_files/'
cmd = 'mkdir ' + TFBS_files_dir
os.system(cmd)

TFBS_file_name = TFBS_bed_file.split('/')[-1]
TFBS_notinTE_file_name = TFBS_file_name.replace('.bed', '_notinTE.bed')
Not_in_TE_TFBS_bed_file = TFBS_files_dir + TFBS_notinTE_file_name
cmd = 'bedtools intersect -v -a '+TFBS_bed_file + ' -b '+TE_TFBS_wb_file + ' > '+ Not_in_TE_TFBS_bed_file
os.system(cmd)

Not_in_TE_TFBS_df = pd.read_csv(Not_in_TE_TFBS_bed_file, sep = '\t', names = ['chromosome', 'geno_start', 'geno_end','TF_name'])
In_TE_TFBS_df = pd.read_csv(TE_TFBS_wb_file, sep = '\t', names = ['chromosome', 'geno_start', 'geno_end','TF_name'])

for bound_TF_name in TF_name_list:
	Not_in_TE_TF_df = Not_in_TE_TFBS_df[Not_in_TE_TFBS_df['TF_name'] == bound_TF_name]
	# print('Not_in_TE_TF_df',Not_in_TE_TF_df)
	Not_in_TF_bed_file = TFBS_files_dir + bound_TF_name + '_notinTE.bed'
	Not_in_TE_TF_df.to_csv(Not_in_TF_bed_file, sep = '\t', header = False, index = False)

	In_TF_df = In_TE_TFBS_df[In_TE_TFBS_df['TF_name'] == bound_TF_name]
	# print('In_TF_df',In_TF_df)
	In_TF_bed_file = TFBS_files_dir + bound_TF_name + '_inTE.bed'
	In_TF_df.to_csv(In_TF_bed_file, sep = '\t', header = False, index = False)

TFBS_coloc_enrich_dir = TFBS_coloc_dir+'Enrichment/'
cmd = 'mkdir ' + TFBS_coloc_enrich_dir
os.system(cmd)


Not_in_TE_count_list = []
Not_in_TE_out_count_list = []
In_TE_count_list = []
In_TE_out_count_list = []
def bedtools_window(TFname):
	i = 0
	while i < len(TF_name_list):
		second_TFname = TF_name_list[i]

		###For the TFBS not in TEs
		Not_in_TE_first_TF_bed = TFBS_files_dir + TFname + '_notinTE.bed'
		Not_in_TE_second_TF_bed = TFBS_files_dir + second_TFname + '_notinTE.bed'

		Not_in_TE_out_temp_file = TFBS_coloc_enrich_dir +'Not_in_TE_out_temp.bed'
		cmd = 'bedtools window -w 500 -u -a '+Not_in_TE_first_TF_bed+' -b '+Not_in_TE_second_TF_bed+' > '+Not_in_TE_out_temp_file
		os.system(cmd)

		Not_in_TE_out_count = Counting_lines_function(Not_in_TE_out_temp_file)
		Not_in_TE_out_count_list.append(int(Not_in_TE_out_count))

		###For the TFBS in TEs
		In_TE_first_TF_bed = TFBS_files_dir + TFname + '_inTE.bed'
		In_TE_second_TF_bed = TFBS_files_dir + second_TFname + '_inTE.bed'

		In_TE_out_temp_file = TFBS_coloc_enrich_dir +'In_TE_out_temp.bed'
		cmd = 'bedtools window -w 500 -u -a '+In_TE_first_TF_bed+' -b '+In_TE_second_TF_bed+' > '+In_TE_out_temp_file
		os.system(cmd)

		In_TE_out_count = Counting_lines_function(In_TE_out_temp_file)
		In_TE_out_count_list.append(int(In_TE_out_count))

		cmd = 'rm '+ In_TE_out_temp_file + ' ' + Not_in_TE_out_temp_file
		os.system(cmd)

		i = i + 1
	return 

j = 0
while j < len(TF_name_list):
	bedtools_window(TF_name_list[j])
	j = j + 1

def Outputting_coloc_count_file(count_list, label_list, raw_output_file):
	count_chunk_list = [count_list[x:x+len(label_list)] for x in range(0, len(count_list), len(label_list))]
	out_df = DataFrame(count_chunk_list,  index = label_list, columns = label_list)
	
	with open(raw_output_file, 'a') as f:
		f.write(out_df.to_csv(header = label_list, index = label_list))
	
	return

NotInTE_dataframe_raw = TFBS_coloc_enrich_dir + TE_name + '_TFBS_NotInTE_colocalisation_count_dataframe_raw.csv'
InTE_dataframe_raw = TFBS_coloc_enrich_dir + TE_name + '_TFBS_InTE_colocalisation_count_dataframe_raw.csv'
Outputting_coloc_count_file(Not_in_TE_out_count_list, TF_name_list, NotInTE_dataframe_raw)
Outputting_coloc_count_file(In_TE_out_count_list, TF_name_list, InTE_dataframe_raw)

###############################################
def Adjusting_for_row(raw_file, label_list, adj_output_file):
	out_df = pd.read_csv(raw_file, header =0, index_col =0)

	out_df_adj = out_df.copy()
	# print(out_df_adj)
	for TF_name in label_list:
		out_df_adj.loc[TF_name] = out_df_adj.loc[TF_name,]/out_df_adj.loc[TF_name, TF_name]
	
	with open(adj_output_file, 'a') as f:
		f.write(out_df_adj.to_csv(header = label_list, index = label_list))

	return

NotInTE_dataframe_adj = TFBS_coloc_enrich_dir + TE_name + '_TFBS_NotInTE_colocalisation_count_dataframe_row_adj.csv'
InTE_dataframe_adj = TFBS_coloc_enrich_dir + TE_name + '_TFBS_InTE_colocalisation_count_dataframe_row_adj.csv'

Adjusting_for_row(NotInTE_dataframe_raw, TF_name_list, NotInTE_dataframe_adj)
Adjusting_for_row(InTE_dataframe_raw, TF_name_list, InTE_dataframe_adj)

NotInTE_df = pd.read_csv(NotInTE_dataframe_raw, header =0, index_col =0)
InTE_df = pd.read_csv(InTE_dataframe_raw, header =0, index_col =0)
NotInTE_df_adj = pd.read_csv(NotInTE_dataframe_adj, header =0, index_col =0)
InTE_df_adj = pd.read_csv(InTE_dataframe_adj, header =0, index_col =0)

fold_change_df = InTE_df_adj.div(NotInTE_df_adj)
fc_adj = TFBS_coloc_enrich_dir + TE_name + '_TFBS_colocalisation_fc.csv'
with open(fc_adj, 'a') as f: 
	f.write(fold_change_df.to_csv(header = TF_name_list, index = TF_name_list))

def twoSampZ(X1, X2, mudiff, n1, n2):
	X1 = float(X1)
	X2 = float(X2)
	n1 = float(n1)
	n2 = float(n2)

	p_hat = (X1+X2)/(n1+n2)
	pooledSE = sqrt(p_hat*(1-p_hat)*(1/n1+1/n2))
	if pooledSE == 0:
		z = 'NA'
		pval = 'NA'
	else:
		z = ((X1/n1 - X2/n2) - mudiff)/pooledSE
		pval = (norm.sf(abs(z)))
	return z, pval

pvalue_df = InTE_df_adj.copy()
zvalue_df = InTE_df_adj.copy()

for first_TF in TF_name_list:
	for second_TF in TF_name_list:
		InTE_obs = InTE_df.loc[first_TF,second_TF]
		InTE_total = InTE_df.loc[first_TF,first_TF]
		NotInTE_obs = NotInTE_df.loc[first_TF,second_TF]
		NotInTE_total = NotInTE_df.loc[first_TF,first_TF]

		z, p = twoSampZ(InTE_obs,NotInTE_obs, 0, InTE_total,NotInTE_total)

		pvalue_df.loc[first_TF,second_TF] = p
		zvalue_df.loc[first_TF,second_TF] = z

pvalue_file = TFBS_coloc_enrich_dir + TE_name + '_TFBS_colocalisation_enrichment_pvalue.csv'
with open(pvalue_file, 'a') as f: 
	f.write(pvalue_df.to_csv(header = TF_name_list, index = TF_name_list))

zvalue_file = TFBS_coloc_enrich_dir + TE_name + '_TFBS_colocalisation_enrichment_zvalue.csv'
with open(zvalue_file, 'a') as f: 
	f.write(zvalue_df.to_csv(header = TF_name_list, index = TF_name_list))

###############################################
###Finding the TSSs in TF-bound TEs
print('Finding the TSSs in TF-bound TEs')
TSS_dir = Work_dir + TEname + '_TFBS_TSS_intersection/'
cmd = 'mkdir ' + TSS_dir
os.system(cmd)

TE_with_TSS_wa = TSS_dir + TEname + '_' + TSS_label + '_TSS_intersection_wa.bed'
cmd = 'bedtools intersect -u -wa -a ' + TE_TFBS_wa_file + ' -b ' +  TSS_bed_file + ' > ' + TE_with_TSS_wa
os.system(cmd)

TE_with_TSS_wb = TSS_dir + TEname + '_' + TSS_label + '_TSS_intersection_wb.bed'
cmd = 'bedtools intersect -u -wa -b ' + TE_TFBS_wa_file + ' -a ' +  TSS_bed_file + ' > ' + TE_with_TSS_wb
os.system(cmd)

TE_with_TSS_wo = TSS_dir + TEname + '_' + TSS_label + '_TSS_intersection_wo.txt'
cmd = 'bedtools intersect -wo -a ' + TE_TFBS_wa_file + ' -b ' +  TSS_bed_file + ' > ' + TE_with_TSS_wo
os.system(cmd)

print('Number of TF-bound TEs with at least one TSS', Counting_lines_function(TE_with_TSS_wa))
print('Number of TSSs that overlap with TF-bound TEs', Counting_lines_function(TE_with_TSS_wb))

print('Done')