#!/usr/bin/python
import pandas as pd
import networkx as nx
import numpy as np
import sys


file_name = sys.argv[1] #"MCF-7_STARR.consolidate.bed.sort.filter.fall_in_region_d.sort.umotif.cre"
window_size = int(sys.argv[2])
step_size = 1

G = nx.Graph()
with open("../Umotif_Interaction_Score.txt") as fin:
    for line in fin:
        line_split = line.strip().split("\t")
        from_node = line_split[0].replace("cluster_", "")
        to_node = line_split[1].replace("cluster_", "")
        weight = float(line_split[2])
        G.add_edge(from_node, to_node, weight=weight)


df_cre = pd.read_csv(file_name, sep="\t", header=None, index_col=None, low_memory=False)
df_cre.iloc[:,-2] = df_cre.iloc[:,-2].astype('str')

df_cre_no_repeat_tfbs = df_cre.iloc[:,:3].join(df_cre.iloc[:, -5:-2])
df_cre_no_repeat_tfbs = df_cre_no_repeat_tfbs.drop_duplicates()
df_cre_no_repeat_tfbs.columns = range(len(df_cre_no_repeat_tfbs.columns))
df_cre_no_repeat_tfbs_ct = df_cre_no_repeat_tfbs.groupby([0,1,2])[3].count().reset_index()
df_cre_no_repeat_tfbs_ct['TFBS_num_in_100bp'] = 100*df_cre_no_repeat_tfbs_ct[3]/(df_cre_no_repeat_tfbs_ct[2] - df_cre_no_repeat_tfbs_ct[1])

df_cre_umotif = df_cre.iloc[:, :3].join(df_cre.iloc[:, -2])
df_cre_umotif.columns = range(len(df_cre_umotif.columns))

def compute_sliding_node_degree(lst, window_size, step_size, G):
    try:
        mean_degree_list = []
        lst = lst.tolist()
        if len(lst) < window_size:
            subnet = G.subgraph(lst)
            # Compute the degree distribution
            degree_list = [degree for _, degree in subnet.degree(weight='weight')]
            #print(degree_list)
            degree_mean = np.mean(degree_list)

        else:
            for i in range(0, len(lst) - window_size + 1, step_size):
                sublist = lst[i:i+window_size]
                subnet = G.subgraph(sublist)
                # Compute the degree distribution
                degree_list = [degree for _, degree in subnet.degree(weight='weight')]
                #print(degree_list)
                degree_mean = np.mean(degree_list)
                #print(degree_mean)
                mean_degree_list.append(degree_mean)
            degree_mean = np.mean(mean_degree_list)
    except:
        degree_mean = 0
            
    return degree_mean




df_cre_umotif_degree = df_cre_umotif.groupby([0,1,2]).agg(
        {3:lambda x:compute_sliding_node_degree(x, window_size=window_size, step_size=1, G=G)}
    ).reset_index()

df_cre_umotif_degree['degree_div_peak_len'] = 100 * df_cre_umotif_degree[3]/(df_cre_umotif_degree[2]-df_cre_umotif_degree[1])

df_cre_umotif_degree_tfbs_ct = df_cre_umotif_degree.merge(df_cre_no_repeat_tfbs_ct, on=[0,1,2])
df_cre_umotif_degree_tfbs_ct = df_cre_umotif_degree_tfbs_ct.rename(mapper={"3_x":"avg_degree", "3_y":"TFBS_num"}, axis=1)

df_cre_umotif_degree_tfbs_ct.to_csv(f"{file_name}.degree_tfbs.{window_size}.csv")
