#!/usr/bin/env python
# coding: utf-8

import networkx as nx    
import argparse
from argparse import ArgumentParser

def get_option(thres_cor, thres_clique_num, thres_cliqueness):
    argparser = ArgumentParser()
    argparser.add_argument('-e', '--edge', type=float,
                           default=thres_cor,
                           help='Correlation coefficient threshold. An edge is drawn if this threshold is exceeded in the gene co-expression network. 0.0 - 1.0 is accepted.')
    argparser.add_argument('-n', '--number', type=int,
                           default=thres_clique_num,
                           help='Minimum number of nodes in clique.')
    argparser.add_argument('-c', '--cliqueness', type=float,
                           default=thres_cliqueness,
                           help='Threshold of cliqueness when merging modules. 0.0 - 1.0 is accepted.')
    return argparser.parse_args()

def preprocessing(df_host_microb, df_host_host, df_microb_microb, thres_cor):
    # Make gene list
    list_host = list((set(df_host_host['gene1']) | set(df_host_host['gene2'])) & (set(df_host_microb['host_gene'])))
    list_microb = list((set(df_microb_microb['gene1']) | set(df_microb_microb['gene2'])) & (set(df_host_microb['microbiome_gene'])))
    
    # Extract gene pairs included in the gene list
    df_host_host = df_host_host[df_host_host['gene1'].isin(list_host)]
    df_host_host = df_host_host[df_host_host['gene2'].isin(list_host)]    
    df_microb_microb = df_microb_microb[df_microb_microb['gene1'].isin(list_microb)]
    df_microb_microb = df_microb_microb[df_microb_microb['gene2'].isin(list_microb)]    
    df_host_microb = df_host_microb[df_host_microb['host_gene'].isin(list_host)]
    df_host_microb = df_host_microb[df_host_microb['microbiome_gene'].isin(list_microb)]
    
    # Extract gene pairs with correlations above a threshold 
    df_host_microb = df_host_microb[abs(df_host_microb['correlation_r'])>thres_cor]
    df_host_host = df_host_host[df_host_host['correlation_r']>thres_cor]
    df_microb_microb = df_microb_microb[df_microb_microb['correlation_r']>thres_cor]
    return df_host_microb, df_host_host, df_microb_microb


def make_graph(df_host_microb, df_host_host, df_microb_microb):
    # Generate tuple of gene pair
    tp_hm = [tuple(x) for x in df_host_microb.iloc[:,[0,1]].values]
    tp_hh = [tuple(x) for x in df_host_host.iloc[:,[0,1]].values]
    tp_mm = [tuple(x) for x in df_microb_microb.iloc[:,[0,1]].values]
    
    # Generate graphs
    G_hm = nx.Graph(tp_hm)
    G_hh = nx.Graph(tp_hh)
    G_mm = nx.Graph(tp_mm)
    G = nx.compose(G_hm,G_hh)
    G = nx.compose(G,G_mm)
    
    # Node list
    list_node_host = list(set(list(set(df_host_microb.loc[:,'host_gene'])) + list(set(df_host_host.loc[:,'gene1'])) + list(set(df_host_host.loc[:,'gene2']))))
    list_node_microbiome = list(set(list(set(df_host_microb.loc[:,'microbiome_gene'])) + list(set(df_microb_microb.loc[:,'gene1'])) + list(set(df_microb_microb.loc[:,'gene2']))))
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
    
    # the number of edges if the community is clique
    num_full_edge = num_nodes_union * (num_nodes_union - 1) / 2
    num_full_edge_host = num_nodes_union_host * (num_nodes_union_host - 1) / 2
    num_full_edge_microb = num_nodes_union_microb * (num_nodes_union_microb - 1) / 2
    num_full_edge_host_microb = num_nodes_union_host * num_nodes_union_microb
    
    # Cliqueness
    ## Host - Host
    if num_full_edge_host == 0: cliqueness_host = 0
    else: cliqueness_host = num_edge_host / num_full_edge_host
    ## Microbiome - Microbiome
    if num_full_edge_microb ==0: cliqueness_microb = 0
    else: cliqueness_microb = num_edge_microb / num_full_edge_microb
    ## Host - Mictobiome
    if num_full_edge_host_microb == 0: cliqueness_host_microbiome = 0
    else: cliqueness_host_microbiome = num_edge_host_microb / num_full_edge_host_microb
    
    return cliqueness_host, cliqueness_microb, cliqueness_host_microbiome


def module_extraction(G, thres_clique_num, thres_clique_rate):
    dict_clique_num = nx.node_clique_number(G)
    list_modules = []
    
    # Sort by number of nodes 
    dict_clique_num = {k: v for k, v in sorted(dict_clique_num.items(), key = lambda x : x[1], reverse=True)}
    
    while len(dict_clique_num)>0:
        if list(dict_clique_num.values())[0]<thres_clique_num: 
            break
        # Find the largest clique
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
                
        # Compute cliqueness between the target clique and all cliques included in the clique list, respectively.
        # If the maximum clique_rate is higher than the threshold, the clique pair is merged. 
        max_clique_rate = 0
        max_clique_num = 0
        for i in range(len(list_modules)):
            cliqueness_host, cliqueness_microb, cliqueness_host_microbiome = compute_cliquness(list_modules[i], list_nodes, G)
            clique_rate_min = min(cliqueness_host, cliqueness_microb, cliqueness_host_microbiome)
            if clique_rate_min > max_clique_rate:
                max_clique_rate = clique_rate_min
                max_clique_num = i
        if max_clique_rate > thres_clique_rate:
            list_modules[max_clique_num] = sorted(list(set(list_modules[max_clique_num]).union(set(list_nodes))))
        else:
            list_modules.append(sorted(list_nodes))
    return list_modules
