#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd 
from sklearn import cluster, preprocessing 
from sklearn import datasets 
import networkx as nx
import networkx.algorithms.community as nx_comm
import copy
import statsmodels.stats.multitest as multi
import argparse
from argparse import ArgumentParser
import sys
from itertools import count

#====================================================================
thres_num_gene_pw_min= 5
thres_num_gene_pw_max = 500
thres_cor = 0.80
thres_clique_num = 20
thres_cliqueness = 0.60
#====================================================================

path_dir_out = sys.argv[1] # Path of output directory
path_spearman_hm = sys.argv[2] # Path of data frame for host-microbiome gene correlation coefficients
path_spearman_hh = sys.argv[3] # Path of data frame for host-host gene correlation coefficients
path_spearman_mm = sys.argv[4] # Path of data frame for microbiome-microbiome gene correlation coefficients
df_hm = pd.read_csv(path_spearman_hm, sep = '\t', header=0, index_col=None, names=['host_gene', 'microbiome_gene', 'correlation_r'])
df_hh = pd.read_csv(path_spearman_hh, sep = '\t', header=0, index_col=None, names=['gene1', 'gene2', 'correlation_r'])
df_mm = pd.read_csv(path_spearman_mm, sep = '\t', header=0, index_col=None, names=['gene1', 'gene2', 'correlation_r'])

# extract gene list including all three dataset
list_host_gene = set(list(df_hh['gene1']) + list(df_hh['gene2']))
list_host_gene = list(set(list_host_gene) & (set(df_hm['host_gene'])))
list_mcrb_gene = set(list(df_mm['gene1']) + list(df_mm['gene2']))
list_mcrb_gene = list((set(list_mcrb_gene)) &  (set(df_hm['microbiome_gene'])))

# extract column with gene including gene list
df_hh = df_hh[df_hh['gene1'].isin(list_host_gene)]
df_hh = df_hh[df_hh['gene2'].isin(list_host_gene)]

df_mm = df_mm[df_mm['gene1'].isin(list_mcrb_gene)]
df_mm = df_mm[df_mm['gene2'].isin(list_mcrb_gene)]

df_hm = df_hm[df_hm['host_gene'].isin(list_host_gene)]
df_hm = df_hm[df_hm['microbiome_gene'].isin(list_mcrb_gene)]

def preprocessing(df_hm, df_hh, df_mm, thres_cor):
    df_hm = df_hm[abs(df_hm['correlation_r'])>thres_cor]
    df_hh = df_hh[df_hh['correlation_r']>thres_cor]
    df_mm = df_mm[df_mm['correlation_r']>thres_cor]

    # Host-Host correlation: extract genes including the correlated host-microbiome gene pair 
    df_hh = df_hh[(df_hh['gene1'].isin(list(set(df_hm.loc[:,'host_gene'])))) & (df_hh['gene2'].isin(list(set(df_hm.loc[:,'host_gene']))))]
    # Microbiome-Microbiome correlation: extract genes including the correlated host-microbiome gene pair 
    df_mm = df_mm[(df_mm['gene1'].isin(list(set(df_hm.loc[:,'microbiome_gene'])))) & (df_mm['gene2'].isin(list(set(df_hm.loc[:,'microbiome_gene']))))]    
    return df_hm, df_hh, df_mm

def make_graph(df_hm, df_hh, df_mm):
    
    # Generate tuple of gene pair
    tp_hm = [tuple(x) for x in df_hm.iloc[:,[0,1]].values]
    tp_hh = [tuple(x) for x in df_hh.iloc[:,[0,1]].values]
    tp_mm = [tuple(x) for x in df_mm.iloc[:,[0,1]].values]
    
    # Generate graphs
    G_hm = nx.Graph(tp_hm)
    G_hh = nx.Graph(tp_hh)
    G_mm = nx.Graph(tp_mm)
    G = nx.compose(G_hm,G_hh)
    G = nx.compose(G,G_mm)
    
    # Node list
    list_node_host = list(set(list(set(df_hm.loc[:,'host_gene'])) + list(set(df_hh.loc[:,'gene1'])) + list(set(df_hh.loc[:,'gene2']))))
    list_node_microbiome = list(set(list(set(df_hm.loc[:,'microbiome_gene'])) + list(set(df_mm.loc[:,'gene1'])) + list(set(df_mm.loc[:,'gene2']))))
    list_nodes = list_node_host + list_node_microbiome
    
    # Attribute list
    list_attr_host = ['Host' for i in  range(len(list_node_host))]
    list_attr_microbiome = ['Microbiome' for i in  range(len(list_node_microbiome))]
    list_attr = list_attr_host + list_attr_microbiome

    # Dictionally 
    group_attr = dict(set(zip(list_nodes,list_attr)))

    # set nodes attributes
    nx.set_node_attributes(G, name="group", values=group_attr)
    
    return G

def compute_cliquness(commun1, commun2, G):
    group_attr = nx.get_node_attributes(G, "group")
    commun1 = frozenset(commun1)
    commun2 = frozenset(commun2)
    # the list of nodes
    commun_union = commun1.union(commun2) # the union set of commun1 and commun2
    commun_union_host = [node for node in commun_union if group_attr[node]=='Host']
    commun_union_microb = [node for node in commun_union if group_attr[node]=='Microbiome']
    # the number of nodes in each list
    num_nodes_union = len(commun_union)
    num_nodes_union_host = len(commun_union_host)
    num_nodes_union_microb = len(commun_union_microb)
    # make subgraph
    G_sub = G.subgraph(commun_union)
    G_sub_host = G.subgraph(commun_union_host)
    G_sub_microb = G.subgraph(commun_union_microb)
    # the list of edges between 'Host' and 'Microbiome' node (bipartite graph)
    edges_host_microb = [[i,j] for i,j in list(G_sub.edges()) if ((group_attr[i]=='Host')&(group_attr[j]=='Microbiome')) |
                         ((group_attr[i]=='Microbiome')&(group_attr[j]=='Host'))]
    # the number of edges
    num_edge = G_sub.number_of_edges()
    num_edge_host = G_sub_host.number_of_edges()
    num_edge_microb = G_sub_microb.number_of_edges()
    num_edge_host_microb = len(edges_host_microb)
    # number of edges if the community is clique
    num_full_edge = num_nodes_union * (num_nodes_union - 1) / 2
    num_full_edge_host = num_nodes_union_host * (num_nodes_union_host - 1) / 2
    num_full_edge_microb = num_nodes_union_microb * (num_nodes_union_microb - 1) / 2
    num_full_edge_host_microb = num_nodes_union_host * num_nodes_union_microb
    
    # Cliqueness
    if num_full_edge_host == 0:cliqueness_host = 0
    else: cliqueness_host = num_edge_host / num_full_edge_host
    if num_full_edge_microb ==0: cliqueness_microb = 0
    else: cliqueness_microb = num_edge_microb / num_full_edge_microb
    if num_full_edge_host_microb == 0: cliqueness_host_microbiome = 0
    else: cliqueness_host_microbiome = num_edge_host_microb / num_full_edge_host_microb
    
    return cliqueness_host, cliqueness_microb, cliqueness_host_microbiome

def module_extraction(G, thres_clique_num, thres_cliqueness):
    dict_clique_num = nx.node_clique_number(G)
    
    list_modules = []
    rec = 0
    while len(dict_clique_num)>0:
        rec+=1
        dict_clique_num = {k: v for k, v in sorted(dict_clique_num.items(), key = lambda x : x[1], reverse=True)}
        if list(dict_clique_num.values())[0]<thres_clique_num: 
            break
        list_candidate = nx.cliques_containing_node(G,list(dict_clique_num.keys())[0])
        for com_num in range(len(list_candidate)):
            list_candidate[com_num] = sorted(list_candidate[com_num])
        list_candidate = list(zip([len(com) for com in list_candidate],list_candidate))
        list_candidate.sort(reverse=True)
        list_candidate = [list_candidate[i][1] for i in range(len(list_candidate))]
        list_nodes = list_candidate[0]
        
        # Delete listed nodes from the dictionary
        for node in list_nodes:
            if node in dict_clique_num: del dict_clique_num[node]
        # Compute clique-rate between the target clique and all clique included in the clique list, respectively.
        # If the maximum clique_rate is higher than the threshold, the clique pair is merged. 
        max_clique_rate = 0
        max_clique_num = 0
        for i in range(len(list_modules)):
            cliqueness_host, cliqueness_microb, cliqueness_host_microbiome = compute_cliquness(list_modules[i], list_nodes, G)
            clique_rate_min = min(cliqueness_host, cliqueness_microb, cliqueness_host_microbiome)
            if clique_rate_min > max_clique_rate:
                max_clique_rate = clique_rate_min
                max_clique_num = i
        if max_clique_rate > thres_cliqueness:
            list_modules[max_clique_num] = sorted(list(set(list_modules[max_clique_num]).union(set(list_nodes))))
        else:
            list_modules.append(sorted(list_nodes))
    return list_modules

# Preprocessing
df_hm_process, df_hh_process, df_mm_process = preprocessing(df_hm, df_hh, df_mm, thres_cor)

# Make Graph
G = make_graph(df_hm_process, df_hh_process, df_mm_process)

# Find modules
list_modules = module_extraction(G, thres_clique_num, thres_cliqueness)

# Output clique member list
path_out = path_dir_out + '/module_list.tsv'    

with open(path_out, 'w', encoding='utf-8', newline='\n') as f:
    for data_slice in list_modules:
        for data in data_slice:
            f.write("%s\t" % data)
        f.write('\n')
