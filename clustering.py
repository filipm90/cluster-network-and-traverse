#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import sys, time
import pandas as pd
import networkx as nx

pd.options.mode.chained_assignment = None #to avoid displaying errors of future obsoleteness in Pandas dataframe

"""
                    *****clustering.py*****

Generation of the network clusters (or components) on the basis of
pairwise node relationships and assignment of cluster number to each
node pair. In a chemical network representation, nodes present compounds
and edges their pairwise relationships (e.g., similarity). Although the
program was developed for the purpose of the chemical space assessment,
the script can be used for any network representation where nodes and edges
(of any origin) are defined.

                    ***Required Packages***

networkx (pip install networkx)
pandas

                        ***How to Run***

python clustering.py test.dat

where:    clustering.py -> script name
          test.dat -> node pair file where first two columns represent node instances forming pairwise relationships per row (e.g., node_1, node_2)

Use file test.dat to run an example clustering!

                        ***Output***

output.dat -> an updated input file.dat with assigned cluster IDs per each edge row
cluster_info.dat -> basic cluster statistics

Developed by Filip Miljković

"""

start = time.time()

file_name = sys.argv[1]

edge_file = pd.read_csv(file_name, sep = "\t")

network_list = list(zip(edge_file["CPD_1"].tolist(), edge_file["CPD_2"].tolist()))

g = nx.Graph()
g.add_edges_from(network_list)

print (nx.info(g))

node_dict = {}
count = 1

for cluster in nx.connected_components(g):
    for node in tuple(cluster):
        node_dict[node] = count
    count += 1

cluster_list = []

for pair in network_list:

    cluster_list.append(node_dict[pair[0]])

edge_file["Cluster_ID"] = cluster_list

edge_file.to_csv("output.dat", sep = "\t", index = False)

print ("Total number of clusters: {}".format(len(set(cluster_list))))

###cluster info

#edges

n_edge = dict(edge_file["Cluster_ID"].value_counts())
n_edge = pd.DataFrame.from_dict(n_edge, orient = "index").reset_index()
n_edge.columns = ["Cluster_ID", "Edges"]

#nodes

n_node = edge_file.groupby("Cluster_ID")["CPD_1", "CPD_2"].agg(set).reset_index()

n_node["Nodes"] = None

for i in range(len(n_node["Cluster_ID"])):
    n_node["Nodes"].iloc[i] = len(n_node[n_node.columns[1]].iloc[i].union(n_node[n_node.columns[2]].iloc[i]))
    
n_node = n_node[["Cluster_ID", "Nodes"]]

#cluster info

cluster_info = n_node.merge(n_edge, on = "Cluster_ID", how = "inner")

cluster_info.to_csv("cluster_info.dat", sep = "\t", index = False)

#time

print("Finished in: {} seconds.".format(round(time.time()-start, 3)))

"""End of the code"""
