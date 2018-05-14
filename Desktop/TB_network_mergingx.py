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
import nxviz
from nxviz.plots import CircosPlot

resistance_folder = "/n/data1/hms/dbmi/farhat/rollingDB/summary_table_resistance2.tsv"
directory =  "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/"
target_folder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/data/0-500/"
filefolder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/files/"
antibiotic = 18 ##Rifampicin-antibiotics_dicts[18]#

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
                
input_file=filefolder+"h37rv_genome_summary.csv"
i=0
names=[]
chrom_start=defaultdict(int)
chrom_end=defaultdict(int)
gene=defaultdict(str)
description=defaultdict(str)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if i >0:
                names.append(row[0])
                gene[row[0]] = row[1]
                chrom_start [row[0]] = row[4]
                chrom_end [row[0]] = row[5]
                description [row[0]] = row[7]
                i+=1
            elif i == 0:
                i+=1



input_file = target_folder+"Single_SNP_prob_500.csv"
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


input_file = target_folder+"r2finals_all_500.csv"
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

input_file = target_folder+"phenotype_SNP_all_500.csv"
phenotype_SNP ={}
r2_description =[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                phenotype_SNP [row[0]] = row[1]
                r2_description.append(row[1])


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
first=0
second=0
value=0
position ={}
node_list=[]
for element in r2_finals:
    if phenotype_SNP [element] == "Resistant":
        if r2_final[element] >= 50:
            split=element.split("_")
            if (split[1]+"_"+split[0]) not in r2maxp:
                if Single_SNP_prob[split[0]] < 0.9 and Single_SNP_prob[split[1]] < 0.9:
                    if int(split[0]) not in network_counts_p:
                        network_counts_p.append(int(split[0]))
                    if int(split[1]) not in network_counts_p:
                        network_counts_p.append(int(split[1]))
                    r2maxp.append(element)
#                    value=r2_final[element]
#                    if split[0] not in resistance_comb:
#                        resistance_comb[split[0]] = value
#                    elif split[0] in resistance_comb:
#                        first=resistance_comb[split[0]]
#                        resistance_comb.update({split[0]: first+value})
#                    if split[1] not in resistance_comb:
#                        resistance_comb[split[1]] = value
#                    elif split[1] in resistance_comb:
#                        second=resistance_comb[split[0]]
#                        resistance_comb.update({split[1]: second+value})
colour_dict={}
group_dict={}
Mb1=[]
Mb2=[]
Mb3=[]
Mb4=[]
Mb5=[]
loc_gene=defaultdict(str)
i=0
for element in range(4411532):
#    if element in (sorted(network_counts_p)):
#        #position[element] = i
#        for ident in names:
#            if int(chrom_start [ident]) <= int(element) <= int(chrom_end[ident]):
#                if gene [ident]:
#                    loc_gene [element] = gene [ident]
#                elif not gene [ident]:
#                    loc_gene [element] = ident
    #node_list.append(tuple((element,(loc_gene[element]),i)))
    if element <= 1000000:
        Mb1.append(element)
        #colour_dict [tuple((element,(loc_gene[element]),i))] = "green"
        #group_dict [tuple((element,(loc_gene[element]),i))] = "0-1Mb"
    elif 1000000 < element <= 2000000:
        Mb2.append(element)
#        colour_dict [tuple((element,(loc_gene[element]),i))] = "black"
#        group_dict [tuple((element,(loc_gene[element]),i))] = "1-2Mb"
    elif 2000000 < element <= 3000000:
        Mb3.append(element)
#        colour_dict [tuple((element,(loc_gene[element]),i))] = "red"
#        group_dict [tuple((element,(loc_gene[element]),i))] = "2-3Mb"   
    elif 3000000 < element <= 4000000:
        Mb4.append(element)
#        colour_dict [tuple((element,(loc_gene[element]),i))] = "blue"
#        group_dict [tuple((element,(loc_gene[element]),i))] = "3-4Mb"
    elif 4000000 < element:
        Mb5.append(element)
#        colour_dict [tuple((element,(loc_gene[element]),i))] = "orange"
#        group_dict [tuple((element,(loc_gene[element]),i))] = "4-4.42"
    #i+=1

G_p.add_nodes_from(Mb1, colour="green", group="0-1Mb")
G_p.add_nodes_from(Mb2, colour="black", group="1-2Mb") 
G_p.add_nodes_from(Mb3, colour="red", group="2-3Mb") 
G_p.add_nodes_from(Mb4, colour="blue", group="3-4Mb") 
G_p.add_nodes_from(Mb5, colour="orange", group="4-4.5Mb")
r2maxpl =[]
r2maxpv=[]<
weight_dict={}
for element in r2_finals:
    if phenotype_SNP [element] == "Resistant":
        if r2_final[element] >= 50:
            split=element.split("_")
            if (split[1]+"_"+split[0]) not in r2maxpl:
                if Single_SNP_prob[split[0]] < 0.9 and Single_SNP_prob[split[1]] < 0.9:
                    G_p.add_edge(int(split[0]), int(split[1]) , weight=r2_final[element])
                    #weight_dict[tuple(((int(split[0]), loc_gene[int(split[0])], position[int(split[0])]), (int(split[1]), loc_gene[int(split[1])], position[int(split[1])]))) = r2_final[element]
                    #edges,weights = zip(*nx.get_edge_attributes(G_p,'weight').items())
                    r2maxpl.append(element)
                    r2maxpv.append(r2_final[element])
                    
#count=Counter(sorted(r2maxpv))
#colors = range(len(count))
#dec={}
#i=0
#for key, value in count.items():
#    dec [key] = i
#    i+=1

color_node=nx.get_node_attributes(G_p,'colour')    
for element in r2maxpl:
    split=element.split("_")
    #G_p [int(split[0])] [int(split[1])] ["color"] = dec[r2_final[element]]
    G_p [int(split[0])] [int(split[1])] ["color"] = color_node [int(split[0])]

fig3 = plt.figure(figsize=(30,30))
c = CircosPlot(G_p, node_color="colour", node_grouping="group", group_order="alphabetically" ,group_label_position="middle", group_label_color=True ,edge_width="weight", edge_cmap=plt.cm.Blues, edge_color="color")
c.draw()
#plt.show()
#fig3.savefig(target_folder+'network_nxviz_colors_single.png', format="PNG", dpi=400)
plt.savefig(target_folder+'circos_500.png', format="PNG", dpi=400)
plt.close()