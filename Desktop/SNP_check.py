#!/usr/bin/env python3

import matplotlib
matplotlib.use("Agg")
import os
import vcf
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pandas import DataFrame
from os import listdir
from os.path import isfile, join
import math
import regex as re
from numpy import array
from collections import defaultdict
from collections import Counter
from sklearn import svm
from mlxtend.plotting import plot_decision_regions
import matplotlib.gridspec as gridspec
import itertools
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from mlxtend.classifier import EnsembleVoteClassifier
from sklearn.naive_bayes import GaussianNB 


resistance_folder = "/n/data1/hms/dbmi/farhat/rollingDB/summary_table_resistance2.tsv"
directory =  "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/"
target_folder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/files/"
filefolder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/"
antibiotic = 18 ##Rifampicin-antibiotics_dicts[18]##

#############################Antibiotics data##################################
phenotypes =[]
isolates =[]
isolate_col =[]
antibiotics =[]
i=0
with open(resistance_folder, "r") as tsvfile:
    for i, row in enumerate(tsvfile):
        while i<1:
            split = row.split("\t")
            antibiotics = split [1:]
            i +=1
            pass
        split= row.split("\t")
        isolate_col.append(split[0])
        i +=1
    isolates = isolate_col [1:]
antibiotics_dicts= [{} for _ in range(len(antibiotics))]
antibiotics_lists=[[] for _ in range(len(antibiotics))]
with open(resistance_folder, "r") as tsvfile:
    for row in tsvfile:
        k=0
        split=row.rstrip("\n").split("\t")
        for k in range (0,(len(split)-1)):                
            if split[k+1] == "R" or split [k+1] == "S":
                antibiotics_dicts[k] [split[0]] = split[k+1]
                antibiotics_lists[k].append(split[0])
                
################################################################################
################for combination increase i<1001#############
file_dict ={}
file_list =[]
name_list =[]
lineages_file={}
lineages_filel=[]
i=0
for root, dir, files in os.walk(directory):
        for name in files:
                if name.endswith('.vcf'):
                    name = name.split(".")
                    names = re.findall(r"(\w+\d+)", name[0])
                    filename = str(names[0])
                    if filename in antibiotics_dicts[antibiotic]:
                        if antibiotics_dicts[antibiotic][filename] in "RS":
                            if i <=500:
                                i+=1
                            elif 500 < i < 1001:
                                file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
                                file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
                                name_list.append(name)
                                lineages_file [name[0]]= directory+name[0]+"/fast-lineage-caller/"+name[0]+".lineage"
                                lineages_filel.append(directory+name[0]+"/fast-lineage-caller/"+name[0]+".lineage")
                                i +=1


input_file = target_folder+"r2finals_all_1000.csv"
r2_final ={}
r2_finals =[]
r2_finalp =[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final [row[0]] = float(row[1])
                r2_finals.append(row[0])
                r2_finalp.append(float(row[1]))

input_file = target_folder+"phenotype_SNP_all_1000.csv"
phenotype_SNP ={}
r2_description =[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                phenotype_SNP [row[0]] = row[1]
                r2_description.append(row[1])
                
input_file = filefolder+"list1/Single_SNP_prob.csv"
Single_SNP_prob ={}
Single_SNP_probs =[]
Single_SNP_probp =[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                Single_SNP_prob [row[0]] = float(row[1])
                Single_SNP_probs.append(row[0])
                Single_SNP_probp.append(float(row[1]))
                

#G_p=nx.Graph()
#G_c=nx.Graph()
#split =[]
#r2maxp=[]
#r2maxc=[]
#single =[]
#network_counts_p =[]
#network_counts_c =[]
#network_counts_p_single =[]
#network_counts_c_single =[]
#resistance_comb={}
#first=0
#second=0
#value=0
##rpoB=[]
#for element in r2_finals:
#    if phenotype_SNP [element] == "Resistant":
#        if r2_final[element] >= 10:
#            split=element.split("_")
##            if "761155" in split:
##                rpoB.append(element)
#            if (split[1]+"_"+split[0]) not in r2maxp:
#                if Single_SNP_prob[split[0]] < 0.9 and Single_SNP_prob[split[1]] < 0.9:
#                    network_counts_p.append(split[0])
#                    network_counts_p.append(split[1])
#                    if split[0] not in network_counts_p:
#                        G_p.add_node(split[0], prob=Single_SNP_prob[split[0]])
#                    if split[1] not in network_counts_p:
#                        G_p.add_node(split[1], prob=Single_SNP_prob[split[1]])
#                    G_p.add_edge(split[0], split[1], fillcolor="red", color="red", weight=r2_final[element])
#                    edges,weights = zip(*nx.get_edge_attributes(G_p,'weight').items())
#                    r2maxp.append(element)
##                    value=r2_final[element]
##                    if split[0] not in resistance_comb:
##                        resistance_comb[split[0]] = value
##                    elif split[0] in resistance_comb:
##                        first=resistance_comb[split[0]]
##                        resistance_comb.update({split[0]: first+value})
##                    if split[1] not in resistance_comb:
##                        resistance_comb[split[1]] = value
##                    elif split[1] in resistance_comb:
##                        second=resistance_comb[split[0]]
##                        resistance_comb.update({split[1]: second+value})
#
#fig3 = plt.figure(figsize=(50,50))            
#nx.draw_circular(G_p,node_color='g', with_labels=True, edge_color=weights, edge_cmap=plt.cm.Blues, node_size=20)
#plt.axis('equal')
#fig3.savefig(target_folder+'network_p_final.pdf', format="PDF", dpi=10000)
#plt.close()


import networkx as nx
import matplotlib.pyplot as plt
G_p=nx.Graph()
G_c=nx.Graph()
split =[]
r2maxp=[]
r2maxc=[]
single =[]
network_counts_p =[]
network_counts_c =[]
network_counts_p_single =[]
network_counts_c_single =[]
resistance_comb={}
results=[]
results2=[]
results3=[]
first=0
second=0
value=0
#rpoB=[]
for element in r2_finals:
    if phenotype_SNP [element] == "Resistant":
        if r2_final[element] >= 10:
            split=element.split("_")
#            if "761155" in split:
#                rpoB.append(element)
            if (split[1]+"_"+split[0]) not in r2maxp:
                if Single_SNP_prob[split[0]] < 0.9 and Single_SNP_prob[split[1]] < 0.9:
                    if int(split[0]) not in network_counts_p:
                        network_counts_p.append(int(split[0]))
                    if int(split[1]) not in network_counts_p:
                        network_counts_p.append(int(split[1]))
                    G_p.add_edge(int(split[0]), int(split[1]), weight=r2_final[element])
                    edges,weights = zip(*nx.get_edge_attributes(G_p,'weight').items())
                    r2maxp.append(element)
#                value=r2_final[element]
#                if split[0] not in resistance_comb:
#                    resistance_comb[split[0]] = value
#                elif split[0] in resistance_comb:
#                    first=resistance_comb[split[0]]
#                    resistance_comb.update({split[0]: first+value})
#                if split[1] not in resistance_comb:
#                    resistance_comb[split[1]] = value
#                elif split[1] in resistance_comb:
#                    second=resistance_comb[split[0]]
#                    resistance_comb.update({split[1]: second+value})

#results = list(map(int, network_counts_p))
#results2= sorted(results)
#results3= list(map(str, results2))
G_p.add_nodes_from(sorted(network_counts_p))
fig3 = plt.figure(figsize=(50,50))            
nx.draw_circular(G_p,node_color='g', with_labels=True, edge_color=weights, edge_cmap=plt.cm.Blues, node_size=20)
plt.axis('equal')
fig3.savefig(target_folder+'network_p_finalsorted_1000.pdf', format="PDF", dpi=10000)
plt.close()

#lineages={}
#lineage=[]
#for file in lineages_filel:
#    try:
#        with open(file, "r") as csvfile:
#            reader = csv.reader(csvfile)
#            for row in reader:
#                name =re.findall(r"\W(\w+)", file.split(".")[0]) [-1]
#                lineages[name] = re.findall("\d", row[0])
#                lineage.append(re.findall('\d', row[0]))
#    #            lineages[name] = re.findall('\d.\d', row[1])
#    #            lineage.append(re.findall('\d.\d', row[1]))
#    except:
#        pass
#    
#s1753519=[]
#s2525722=[]
#s2881597=[]
#s3415180=[]
#s131174=[]
#
#
#
#names=[]
#for key,val in file_dict.items():
#    if key in lineages:
#        pos = []
#        with open(val) as vcffile:
#            pos=[record.POS for record in vcf.Reader(vcffile) if record.is_snp]
#        if 1753519 in pos:
#            s1753519.append(lineages[key])
#        elif 2525722 in pos:
#            s2525722.append(lineages[key])
#        elif 2881597 in pos:
#            s2881597.append(lineages[key])
#        elif 3415180 in pos:
#            s3415180.append(lineages[key])
#        elif 131174 in pos:
#            s131174.append(ineages[key])
#        
#d1753519= Counter(s1753519)
#d2525722= Counter(s2525722)
#d2881597= Counter(s2881597)
#d3415180= Counter(s3415180)
#d131174= Counter(s131174)
#
#with open (target_folder+"d1753519.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in d1753519.items():
#        writer.writerow([key, value])    
#
#with open (target_folder+"d2525722.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in d2525722.items():
#        writer.writerow([key, value])    
#        
#with open (target_folder+"d2881597.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in d2881597.items():
#        writer.writerow([key, value])    
#        
#with open (target_folder+"d3415180.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in d3415180.items():
#        writer.writerow([key, value])    
#        
#with open (target_folder+"d131174.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in d131174.items():
#        writer.writerow([key, value])    