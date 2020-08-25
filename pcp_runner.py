#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import os, sys, time
import networkx as nx
import pandas as pd
from itertools import combinations

pd.options.mode.chained_assignment = None #to avoid displaying errors of future obsoleteness in Pandas dataframe

"""
                    *****pcp_runner.py*****

Program asks user for input cluster X which will be subjected to pathway analysis
(see cluster_info.dat generated following clustering.py for further cluster information).
Node promiscuity (node_promiscuity) argument is used to define the minimal promiscuity degree (PD)
value that will be used to define source and target nodes for pathways. Pathways are identified,
extracted, and ranked according to predefined criteria. Ranked pathways of cluster X are saved in
./X/pathways_X.dat where X is the cluster number corresponding to information in cluster_info.dat.

                    ***Required Packages***

networkx (pip install networkx)
pandas

                        ***How to Run***

python pcp_runner.py output.dat cluster_info.dat node_promiscuity

where:    pcp_runner.py -> script name
          output.dat -> output file following clustering.py
          cluster_info.dat -> cluster info file following clustering.py
          node_promiscuity -> int number denoting PD of nodes used as source and target in pathway definition

                        ***Output***

./X/cluster_X.dat -> a slice of output.dat with compound pairs coming from a cluster X
./X/pathways_X.dat -> ranked pathways of chosen cluster X with additional pathway information

Developed by Filip Miljković

"""

def mkdir(folder):
    
    """
    Making a new directory for cluster analysis;
    """
    
    new_folder = os.path.join(os.getcwd(), folder)
    
    try:
        
        if not os.path.exists(new_folder):
            os.makedirs(new_folder)
    
    except OSError:
        print ('Error: Creating directory. ' + folder)
        
    return new_folder

def getting_max_uni_prom(d):
    
    """
    Extraction of a shortest path for a node pair with maximum cumulative deltaPD; for several such exist paths, a single instance with same set of promiscuous nodes is retained;
    """
    
    maxv = max(d.values())
    keep_max = {k:v for (k,v) in d.items() if v == maxv}
    uni_prom = {k[0::2]:k for (k,v) in keep_max.items()}
    return {k:v for (k,v) in keep_max.items() if k in uni_prom.values()}

def pairwise_consec(coll):
    
    """
    Edge extraction per pathway;
    """
    
    return zip(coll,coll[1:])

def sort_within(list_of_tuples):
    
    """
    Tuple sorting;
    """
    
    return set(tuple(sorted(t)) for t in list_of_tuples)

"""Start of the code"""

#Checking availability of a Cluster ID and choosing the cluster for further analysis

start = time.time()

pc_file = sys.argv[1]
stat = sys.argv[2]
node_choice = int(sys.argv[3]) #Choose minimal PD of your promiscuous node that will be used to form pathways (e.g., in publication it was set to 6)

stat = pd.read_csv(stat, sep = "\t") #loading cluster statistics file

while True: #checking the availability of the cluster    
    
    q = int(input("Which cluster ID would you like to analyze?\n"))
    qList = stat["Cluster_ID"].tolist()
    
    if q in qList:
        row = int(stat.loc[stat["Cluster_ID"] == q].index.values)
        print ("You chose Cluster ID {} with {} compounds and {} PCs.".format(stat["Cluster_ID"].iloc[row], stat["Nodes"].iloc[row], stat["Edges"].iloc[row]))
        break
        
    else:
        print ("The Cluster ID doesn't exist. You can check available clusters in the file named test_cluster_info.dat. Please try another one...")

#Filtering the cluster file and saving it separately for future use

path = mkdir(str(q))
pc_to_save = "cluster_{}.dat".format(q)
total_path = os.path.join(path, pc_to_save)
pc = pd.read_csv(pc_file, sep = "\t")
pc_cluster = pc[pc["Cluster_ID"] == q]
pc_cluster.to_csv(total_path, sep = "\t", index = False)

pc_cluster["delta_PD"] = abs(pc_cluster["PD_1"] - pc_cluster["PD_2"]) #calculating delta_PD value

#Constructing a graph from the edges

print ("Constructing a network...")

F = nx.Graph(cid = q)

pc_dict = list(zip(pc_cluster["CPD_1"].tolist(), pc_cluster["CPD_2"].tolist(), pc_cluster["delta_PD"].tolist()))
pc_dict = {edge[0:2]: edge[-1] for edge in pc_dict}

inh_dict = list(zip(pc_cluster["CPD_1"], pc_cluster["PD_1"])) + list(zip(pc_cluster["CPD_2"], pc_cluster["PD_2"]))
inh_dict = {node[0]: node[1] for node in inh_dict}

F.add_edges_from(pc_dict.keys()) #forming a network from a chosen cluster

del pc_cluster

print ("Making all pairwise combinations of non-terminal nodes with high promiscuity and calculating the shortest paths...")
    
inner_prom_nodes = filter(lambda node: (len(nx.edges(F, node)) >= 2 and inh_dict[node] >= node_choice), F.nodes()) #only promiscuous compounds with hub-forming abilities

#Constructing all shortest paths for node pairwise combinations

path_set = set()

for i in combinations(inner_prom_nodes, 2):
    
    path_list = []    
    path_list = [tuple(path_tup) for path_tup in nx.all_shortest_paths(F, source = i[0], target = i[1])]
    path_set.add(tuple(path_list))

#Pathway caluculation; saving the output;

print ("Calculating the pathways...")

obj_path = "pathways_{}.dat".format(q)
path_to_obj = os.path.join(path, obj_path)

col_list = ["Pathway No", "Pathway", "Nodes", "No of Edges", "tot_deltaPD"]
output = pd.DataFrame(columns=col_list)

no = 0

pc_dict_sorted = {tuple(sorted(k)):v for (k,v) in pc_dict.items()}

for pair_paths in path_set:
    
    pathway_dictionary = {}
    
    for p_path in pair_paths:
        
        sum_dpd = 0; final = []; temp_dict = {} #defining variables
        pd_list = pairwise_consec(p_path)
        pd_list = [tuple(sorted(tup)) for tup in pd_list]
        temp_dict = {k: v for (k, v) in pc_dict_sorted.items() if k in pd_list}
        sum_dpd = sum(temp_dict.values())
        pathway_dictionary[p_path] = sum_dpd
        
    for k, v in getting_max_uni_prom(pathway_dictionary).items():
        no += 1      
        
        temp_df = pd.DataFrame([[no, k, len(k), len(set(nx.edges(F,k[0::2]))), v]], columns=col_list)
        output = output.append(temp_df, ignore_index=True)

print ("Eliminating redundant pathways...")

#New part of the code, eliminating redundant pathways and saving new file on that way

path_dict = {k: set(v) for k, v in zip(output["Pathway No"], output["Pathway"])} #making dictionary with pathway no as keys and pathway as values

key_list = [k for k in sorted(path_dict, key = lambda k: len(path_dict[k]), reverse = True)]

redundant = set(i[1] for i in combinations(key_list, 2) if (len(path_dict[i[0]]) > len(path_dict[i[1]])) and (path_dict[i[1]].issubset(path_dict[i[0]])))

output = output[~output["Pathway No"].isin(redundant)]

del key_list, redundant

#calculating the tuple of average ranking, sorting the values within the tuples and then sorting the table according to tuple sorting

output["Nodes Rank"] = output["Nodes"].rank(ascending = False)
output["No of Edges Rank"] = output["No of Edges"].rank(ascending = False)
output["tot_deltaPD Rank"] = output["tot_deltaPD"].rank(ascending = False)

output["Fusion Rank"] = list(zip(output["Nodes Rank"], output["No of Edges Rank"], output["tot_deltaPD Rank"]))
output["Fusion Rank"] = [tuple(sorted(x)) for x in output["Fusion Rank"]]
output = output.sort_values(["Fusion Rank"], ascending = True)

output = output[["Pathway No", "Pathway", "Nodes", "No of Edges", "tot_deltaPD", "Fusion Rank"]] #retaining valid columns

output.to_csv(path_to_obj, sep = "\t", index = False) #overwriting previous file with a save

#time

print("Finished in: {} seconds.".format(round(time.time()-start, 3)))

"""End of the code"""