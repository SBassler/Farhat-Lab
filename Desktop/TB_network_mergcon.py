#!/usr/bin/env python3

import os
import vcf
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pandas import DataFrame
import math
import regex as re
from collections import Counter

resistance_folder = "/n/data1/hms/dbmi/farhat/rollingDB/summary_table_resistance2.tsv"
directory =  "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/"
target_folder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/"
filefolder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/"
antibiotic = 18 ##Rifampicin-antibiotics_dicts[18]##

##############################Antibiotics data##################################
#phenotypes =[]
#isolates =[]
#isolate_col =[]
#antibiotics =[]
#i=0
#with open(resistance_folder, "r") as tsvfile:
#    for i, row in enumerate(tsvfile):
#        while i<1:
#            split = row.split("\t")
#            antibiotics = split [1:]
#            i +=1
#            pass
#        split= row.split("\t")
#        isolate_col.append(split[0])
#        i +=1
#    isolates = isolate_col [1:]
#antibiotics_dicts= [{} for _ in range(len(antibiotics))]
#antibiotics_lists=[[] for _ in range(len(antibiotics))]
#with open(resistance_folder, "r") as tsvfile:
#    for row in tsvfile:
#        k=0
#        split=row.rstrip("\n").split("\t")
#        for k in range (0,(len(split)-1)):                
#            if split[k+1] == "R" or split [k+1] == "S":
#                antibiotics_dicts[k] [split[0]] = split[k+1]
#                antibiotics_lists[k].append(split[0])
#                
################################################################################
#
#file_dict ={}
#file_list =[]
#name_list =[]
#i=0
#for root, dir, files in os.walk(directory):
#        for name in files:
#                if name.endswith('.vcf'):
#                    name = name.split(".")
#                    names = re.findall(r"(\w+\d+)", name[0])
#                    filename = str(names[0])
#                    if filename in antibiotics_dicts[antibiotic]:
#                        if antibiotics_dicts[antibiotic][filename] == "R" or "S":
#                            if i < 501:
#                                file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
#                                file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
#                                name_list.append(name)
#                                i +=1
                                

        
input_file = filefolder+"Single_SNP_prob.csv"
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

input_file = filefolder+"Multi_SNP_prob.csv"
Multi_SNP_prob ={}
Multi_SNP_probs =[]
Multi_SNP_probp =[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                if float(row[1]) > 0:
                    Multi_SNP_prob [row[0]] = float(row[1])
                    Multi_SNP_probs.append(row[0])
                    Multi_SNP_probp.append(float(row[1]))


                
###############################################################################
filefolder = filefolder +"files/r2/"
allfiles = [f for f in listdir(filefolder) if isfile(join(filefolder, f))]
os.chdir(filefolder)
r2_final={}
r2_finals =[]
r2_finalp =[]
for file in allfiles:
    with open(input_file, "r") as csvfile:
            reader = csv.reader(csvfile, delimiter = "\t")
            for row in reader:
                if row:
                    r2_final [row[0]] = float(row[1])
                    r2_finals.append(row[0])
                    r2_finalp.append(float(row[1]))
                
##########################Phenotype classification##############################
#filefolder = filefolder +"files/description/"
#allfiles = [f for f in listdir(filefolder) if isfile(join(filefolder, f))]
#os.chdir(filefolder)
#description={}
#for file in allfiles:
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    description [row[0]] = row[1]
#
#r2_description =[]
#phenotype_SNP={}
#for element in r2_finals:
#    if description[element] in ["a","b", "c", "d"]:
#        phenotype_SNP [element] = "Resistant"
#        r2_description.append("Resistant")  
#    elif description [element] in ["e","f", "g", "h"]:
#        phenotype_SNP [element] = "Sensitive"
#        r2_description.append("Sensitive")  
#        
#with open (target_folder+"phenotype_SNP.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in phenotype_SNP.items():
#        writer.writerow([key, value]) 
        
filefolder = filefolder +"files/description/"
allfiles = [f for f in listdir(filefolder) if isfile(join(filefolder, f))]
os.chdir(filefolder)        
phenotype_SNP ={}
r2_description =[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                phenotype_SNP [row[0]] = row[1]
                r2_description.append(row[1])        
############################Manhatten p values#################################    
df = DataFrame({"SNP":r2_finals, "pvalue":r2_finalp, "phenotype": r2_description})
df.phenotype = df.phenotype.astype('category')
df.phenotype = df.phenotype.cat.set_categories(['Resistant', "Sensitive"], ordered=True)
df = df.sort_values('phenotype')

# How to plot gene vs. -log10(pvalue) and colour it by chromosome?
df['minuslog10pvalue'] = -np.log10(df.pvalue)
df['ind'] = range(len(df))
df_grouped = df.groupby(('phenotype'))

fig = plt.figure()
ax = fig.add_subplot(111)
colors = ['red','green']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1   ] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
ax.set_ylim([0, 5])
ax.set_xlabel('phenotype')
fig.savefig(target_folder+'Manhatten.png')

###############SNPs highly connected in pathogenic network#####################
#r2_finalp.sort()
G_p=nx.Graph()
G_c=nx.Graph()
split =[]
single =[]
network_counts_p =[]
network_counts_c =[]
network_counts_p_single =[]
network_counts_c_single =[]
for element in r2_finals:
    if phenotype_SNP [element] == "Resistant":
        if -np.log10(r2_final[element]) > 3:
            split=element.split("_")
            network_counts_p.append(split[0])
            network_counts_p.append(split[1])
            G_p.add_node(split[0], fillcolor="red", color="red", prob=Single_SNP_prob[split[0]])
            G_p.add_node(split[1], fillcolor="red", color="red", prob=Single_SNP_prob[split[1]])
            G_p.add_edge(split[0], split[1], fillcolor="red", color="red", weight=r2_final[element])
        
    elif phenotype_SNP [element] == "Sensitive":
        if -np.log10(r2_final[element]) > 3:
            split=element.split("_")
            network_counts_c.append(split[0])
            network_counts_c.append(split[1])
            G_c.add_node(split[0], fillcolor="green", color="green", prob=Single_SNP_prob[split[0]])
            G_c.add_node(split[1], fillcolor="green", color="green", prob=Single_SNP_prob[split[1]])
            G_c.add_edge(split[0], split[1], fillcolor="green", color="green", weight=r2_final[element])

network_counts_p.sort()
for element in network_counts_p:
    if element not in network_counts_p_single:
        network_counts_p_single.append(element)

network_counts_c.sort()
for element in network_counts_p:
    if element not in network_counts_c_single:
        network_counts_c_single.append(element)

nx.draw_networkx(G_p, with_labels=True, with_fillcolors=True, with_colors=True)
G_p.savefig(target_folder+'G_p network.png')
nx.draw_networkx(G_c, with_labels=True, with_fillcolors=True, with_colors=True)
G_c.savefig(target_folder+'G_c network.png')