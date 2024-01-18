#!/usr/bin/env python
# coding: utf-8

import pandas as pd 
import meign_modules
import sys

#====================================================================
# default parameters
thres_cor = 0.80
thres_clique_num = 20
thres_cliqueness = 0.60
#====================================================================

if __name__ == '__main__':
    args = get_option(thres_cor, thres_clique_num, thres_cliqueness)

path_dir_out = sys.argv[1] # Path of output directory
path_hm = sys.argv[2] # Path of data frame for host-microbiome gene correlation coefficients
path_hh = sys.argv[3] # Path of data frame for host-host gene correlation coefficients
path_mm = sys.argv[4] # Path of data frame for microbiome-microbiome gene correlation coefficients
df_host_microb = pd.read_csv(path_spearman_hm, sep = '\t', header=0, index_col=None, names=['host_gene', 'microbiome_gene', 'correlation_r'])
df_host_host = pd.read_csv(path_spearman_hh, sep = '\t', header=0, index_col=None, names=['gene1', 'gene2', 'correlation_r'])
df_microb_microb = pd.read_csv(path_spearman_mm, sep = '\t', header=0, index_col=None, names=['gene1', 'gene2', 'correlation_r'])

# Preprocessing
df_host_microb, df_host_microb, df_microb_microb = meign_modules.preprocessing(df_host_microb, df_host_host, df_microb_microb, thres_cor)

# Make Graph
G = meign_modules.make_graph(df_host_microb, df_host_microb, df_microb_microb)

# Find modules
list_modules = meign_modules.module_extraction(G, thres_clique_num, thres_cliqueness)

# Output clique member list
path_out = path_dir_out + '/module_list.tsv'    

with open(path_out, 'w', encoding='utf-8', newline='\n') as f:
    for data_slice in list_modules:
        for data in data_slice:
            f.write("%s\t" % data)
        f.write('\n')
