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
antibiotic = 18 ##rifampicin-antibiotics_dicts[18]##

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
            if split[k+1] in "RS":
                antibiotics_dicts[k] [split[0]] = split[k+1]
                antibiotics_lists[k].append(split[0])
                
###############################################################################
file_dict ={}
file_list =[]
name_list =[]
i=0
for root, dir, files in os.walk(directory):
        for name in files:
                if name.endswith('.vcf'):
                    name = name.split(".")
                    names = re.findall(r"(\w+\d+)", name[0])
                    filename = str(names[0])
                    if filename in antibiotics_dicts[antibiotic]:
                        if antibiotics_dicts[antibiotic][filename] in "RS":
                            if i < 501:
                                file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
                                file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
                                name_list.append(name)
                                i +=1
                                
#########Single_SNP_probability and Multi_SNP_probability (1 hour)############
pos=[]
snp1 =[]
snp2=[]
allcombination=[]
allpos=[]
for file in file_list:
    pos=[]
    with open(file) as vcffile:
        vcfReader = vcf.Reader(vcffile)
        for record in vcfReader:
            pos.append(record.POS)
            allpos.append(record.POS)
        i=0
        for first in pos[:-1]:
            snp1=first
            for second in pos[i+1:]:
                snp2=second
                allcombination.append(str(first)+"_"+str(second))
            i+=1
            

Single_SNP_prob ={}
Single_SNP_proba ={}
Single_SNP_probs =[]
Single_SNP_probp =[]
counts = Counter(allpos)
Single_SNP_proba=dict(counts)
for element in Single_SNP_proba:
    Single_SNP_prob[str(element)] = (float(Single_SNP_proba[element]/len(file_list)))
    Single_SNP_probs.append(str(element))
    Single_SNP_probp.append(float(Single_SNP_proba[element]/len(file_list)))
    
with open (target_folder+"Single_SNP_prob.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in Single_SNP_prob.items():
        writer.writerow([key, value])

Multi_SNP_prob ={}
Multi_SNP_proba ={}
Multi_SNP_probs =[]
Multi_SNP_probp =[]
counts = Counter(allcombination)
Multi_SNP_proba=dict(counts)
for element in Multi_SNP_proba:
    if (float(Multi_SNP_proba[element])/len(file_list)) >0:
        Multi_SNP_prob[str(element)] = (float(Multi_SNP_proba[element])/len(file_list))
        Multi_SNP_probs.append(str(element))
        Multi_SNP_probp.append((float(Multi_SNP_proba[element])/len(file_list)))
    
with open (target_folder+"Multi_SNP_prob.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in Multi_SNP_prob.items():
        writer.writerow([key, value])
        


LD = {}
first = []
for key, value in Multi_SNP_prob.items():
    first = re.findall(r"(\d+)",key)
    LD [key] = value-(Single_SNP_prob[first[0]]*Single_SNP_prob[first[1]])
    
with open (target_folder+"LD.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in LD.items():
        writer.writerow([key, value])

pearson= {}
for key, value in LD.items():
    first = re.findall(r"(\d+)",key)
    if (Single_SNP_prob[first[0]]*(1-Single_SNP_prob[first[0]])*Single_SNP_prob[first[1]]*(1-Single_SNP_prob[first[1]]))>0:
        pearson [key] = value/math.sqrt(Single_SNP_prob[first[0]]*(1-Single_SNP_prob[first[0]])*Single_SNP_prob[first[1]]*(1-Single_SNP_prob[first[1]]))
        
with open (target_folder+"pearson.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in pearson.items():
        writer.writerow([key, value])    
    
pearson2={}
for key, value in pearson.items():
    pearson2 [key] = value**2
    
with open (target_folder+"pearson2.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in pearson2.items():
        writer.writerow([key, value])
