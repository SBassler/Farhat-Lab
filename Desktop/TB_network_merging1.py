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
target_folder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/list1/"
filefolder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/"
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
            if split[k+1] == "R" or split [k+1] == "S":
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
                        if antibiotics_dicts[antibiotic][filename] == "R" or "S":
                            if i < 501:
                                file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
                                file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
                                name_list.append(name)
                                i +=1
                                
##########Single_SNP_probability and Multi_SNP_probability (1 hour)############
#pos=[]
#snp1 =[]
#snp2=[]
#allcombination=[]
#allpos=[]
#for file in file_list:
#    pos=[]
#    with open(file) as vcffile:
#        vcfReader = vcf.Reader(vcffile)
#        for record in vcfReader:
#            pos.append(record.POS)
#            allpos.append(record.POS)
#        i=0
#        for first in pos[:-1]:
#            snp1=first
#            for second in pos[i+1:]:
#                snp2=second
#                allcombination.append(str(first)+"_"+str(second))
#            i+=1
#            
#
#Single_SNP_prob ={}
#Single_SNP_proba ={}
#Single_SNP_probs =[]
#Single_SNP_probp =[]
#counts = Counter(allpos)
#Single_SNP_proba=dict(counts)
#for element in Single_SNP_proba:
#    Single_SNP_prob[element] = float((Single_SNP_proba[element]/len(file_list)))
#    Single_SNP_probs.append(str(element))
#    Single_SNP_probp.append(float((Single_SNP_proba[element]/len(file_list))))
#    
#with open (target_folder+"Single_SNP_prob.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in Single_SNP_prob.items():
#        writer.writerow([key, value])
#
#Multi_SNP_prob ={}
#Multi_SNP_proba ={}
#Multi_SNP_probs =[]
#Multi_SNP_probp =[]
#counts = Counter(allcombination)
#Multi_SNP_proba=dict(counts)
#for element in Multi_SNP_proba:
#    Multi_SNP_prob[element] = float((Multi_SNP_proba[element]/len(file_list)))
#    Multi_SNP_probs.append(str(element))
#    Multi_SNP_probp.append(float((Multi_SNP_proba[element]/len(file_list))))
#    
#with open (target_folder+"Multi_SNP_prob.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in Multi_SNP_prob.items():
#        writer.writerow([key, value])
        
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

#LD = {}
#first = []
#for key, value in Multi_SNP_prob.items():
#    first = re.findall(r"(\d+)",key)
#    LD [key] = value-(Single_SNP_prob[first[0]]*Single_SNP_prob[first[1]])
#    
#with open (target_folder+"LD.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in LD.items():
#        writer.writerow([key, value])
#
#pearson= {}
#for key, value in LD.items():
#    first = re.findall(r"(\d+)",key)
#    if (Single_SNP_prob[first[0]]*(1-Single_SNP_prob[first[0]])*Single_SNP_prob[first[1]]*(1-Single_SNP_prob[first[1]]))>0:
#        pearson [key] = value/math.sqrt(Single_SNP_prob[first[0]]*(1-Single_SNP_prob[first[0]])*Single_SNP_prob[first[1]]*(1-Single_SNP_prob[first[1]]))
#        
#with open (target_folder+"pearson.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in pearson.items():
#        writer.writerow([key, value])    
#    
#pearson2={}
#for key, value in pearson.items():
#    pearson2 [key] = value**2
#    
#with open (target_folder+"pearson2.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in pearson2.items():
#        writer.writerow([key, value])

                
##################Phenotype correlations for different drugs###################
clone_key ={}
count_list_a = {}
count_list_b= {}
count_list_c= {}
count_list_d= {}
count_list_e= {}
count_list_f= {}
count_list_g = {}
count_list_h= {}

names=[]
i=0
resistant = 0
sensitive = 0
clone_lists = [[] for _ in range(len(file_dict))]
for key,val in file_dict.items():
    names = re.findall(r"(\w+\d+)", key)
    filename = str(names[0])
    if filename in antibiotics_dicts[antibiotic]:
        if antibiotics_dicts[antibiotic][filename] == "R" or "S":
            if antibiotics_dicts[antibiotic][filename] == "R":
                resistant +=1
            elif antibiotics_dicts[antibiotic][filename] == "S":
                sensitive +=1
            pos = []
            with open(val) as vcffile:
                vcfReader = vcf.Reader(vcffile)
                for record in vcfReader:
                    if record.is_snp:   
                        pos.append(record.POS)
                clone_lists[i]=pos
                clone_key [i] = antibiotics_dicts[antibiotic][filename]
                i+=1  
                
                
clone_listsc = [x for x in clone_lists if x]
count=0
resistance_prob = (resistant/sensitive)
counts_isolates = (len(clone_listsc)) 
length = math.ceil((len(Multi_SNP_probs))/10)
for element in Multi_SNP_probs[:length]:
    if count == (0.5*length):
        with open (target_folder+"count_list_a.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_a.items():
                writer.writerow([key, value])
   
        with open (target_folder+"count_list_b.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_b.items():
                writer.writerow([key, value])  
        
        with open (target_folder+"count_list_c.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_c.items():
                writer.writerow([key, value])

        with open (target_folder+"count_list_d.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_d.items():
                writer.writerow([key, value])

        with open (target_folder+"count_list_e.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_e.items():
                writer.writerow([key, value])
       
        with open (target_folder+"count_list_f.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_f.items():
                writer.writerow([key, value])    
    
        with open (target_folder+"count_list_g.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_g.items():
                writer.writerow([key, value])
       
        with open (target_folder+"count_list_h.csv", "w") as csv_file:
            writer=csv.writer(csv_file, delimiter ="\t")
            for key, value in count_list_h.items():
                writer.writerow([key, value])        
    pos = []
    a=0
    b=0
    c=0
    d=0
    e=0
    f=0
    g=0
    h=0
    i=0
    count_list ={}
    split=element.split("_")
    for clone in clone_listsc:
        if clone_key [i] =="R":
            if int(split[0]) in clone:
                if int(split[1]) in clone:
                    a+=1
                elif int(split[1]) not in clone:
                    b+=1
            elif int(split[0]) not in clone:
                if int(split[1]) in clone:
                    c+=1
                elif int(split[1]) not in clone:
                    d+=1
        elif clone_key [i] == "S":
            if int(split[0]) in clone:
                if int(split[1]) in clone:
                    e+=1
                elif int(split[1]) not in clone:
                    f+=1
            elif int(split[0]) not in clone:
                if int(split[1]) in clone:
                    g+=1
                elif int(split[1]) not in clone:
                    h+=1
        i+=1           
    count_list_a [element] = a/counts_isolates
    count_list_b [element] = b/counts_isolates
    count_list_c [element] = c/counts_isolates
    count_list_d [element] = d/counts_isolates
    count_list_e [element] = e/counts_isolates
    count_list_f [element] = f/counts_isolates
    count_list_g [element] = g/counts_isolates
    count_list_h [element] = h/counts_isolates
    count +=1
    
with open (target_folder+"count_list_a2.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in count_list_a.items():
        writer.writerow([key, value])
   
with open (target_folder+"count_list_b2.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in count_list_b.items():
        writer.writerow([key, value])  
        
with open (target_folder+"count_list_c2.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in count_list_c.items():
        writer.writerow([key, value])

with open (target_folder+"count_list_d2.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in count_list_d.items():
       writer.writerow([key, value])

with open (target_folder+"count_list_e2.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in count_list_e.items():
       writer.writerow([key, value])
       
with open (target_folder+"count_list_f2.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in count_list_f.items():
       writer.writerow([key, value])    
    
with open (target_folder+"count_list_g2.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in count_list_g.items():
       writer.writerow([key, value])
       
with open (target_folder+"count_list_h2.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in count_list_h.items():
       writer.writerow([key, value])    
    
#############calculate R for SNP pairs and phenotype###########################
r_final={}
r_pair_final ={}
r_pair_single_final ={}
r_finallist =[]
r_finalp =[]
r_pair_finalp =[]
r_pair_single_finalp =[]
description ={}
description_pair ={}
description_pair_single ={}
a_d =0
b_d=0
c_d=0
d_d=0
e_d=0
f_d=0
g_d=0
h_d=0
a_dv=0
b_dv=0
c_dv=0
d_dv=0
e_dv=0
f_dv=0
g_dv=0
h_dv=0
split =[]
for element in Multi_SNP_probs[:length]:
    split=element.split("_")
    counts = []
    counts_pair = []
    counts_pair_single = []
    keys ={}
    keys_pair={}
    keys_pair_single={}
    a = count_list_a[element]-(Multi_SNP_prob[element]*(resistance_prob))
    b = count_list_b[element]-(Multi_SNP_prob[element]*(resistance_prob))
    c = count_list_c[element]-(Multi_SNP_prob[element]*(resistance_prob))
    d = count_list_d[element]-(Multi_SNP_prob[element]*(resistance_prob))
    e = count_list_e[element]-(Multi_SNP_prob[element]*(1-resistance_prob))
    f = count_list_f[element]-(Multi_SNP_prob[element]*(1-resistance_prob))
    g = count_list_g[element]-(Multi_SNP_prob[element]*(1-resistance_prob))
    h = count_list_h[element]-(Multi_SNP_prob[element]*(1-resistance_prob))
    a_dv=(float(Multi_SNP_prob[element]))*(1-(float(Multi_SNP_prob[element])))*(resistance_prob)*(1-resistance_prob)
    b_dv=(float(Multi_SNP_prob[element]))*(1-(float(Multi_SNP_prob[element])))*(resistance_prob)*(1-resistance_prob)
    c_dv=(float(Single_SNP_prob[split[0]])*(1-(float(Single_SNP_prob[split[1]]))))*(1-(float(Single_SNP_prob[split[0]])*(1-(float(Single_SNP_prob[split[1]]))))*(resistance_prob)*(1-resistance_prob))
    d_dv=(float(Single_SNP_prob[split[0]])*(1-(float(Single_SNP_prob[split[1]]))))*(1-(float(Single_SNP_prob[split[0]])*(1-(float(Single_SNP_prob[split[1]]))))*(resistance_prob)*(1-resistance_prob))
    e_dv=(float(Single_SNP_prob[split[1]])*(1-(float(Single_SNP_prob[split[0]]))))*(1-(float(Single_SNP_prob[split[1]])*(1-(float(Single_SNP_prob[split[0]]))))*(resistance_prob)*(1-resistance_prob))
    f_dv=(float(Single_SNP_prob[split[1]])*(1-(float(Single_SNP_prob[split[0]]))))*(1-(float(Single_SNP_prob[split[1]])*(1-(float(Single_SNP_prob[split[0]]))))*(resistance_prob)*(1-resistance_prob))
    g_dv=((1-(float(Single_SNP_prob[split[0]])))*(1-(float(Single_SNP_prob[split[1]]))))*(1-((1-(float(Single_SNP_prob[split[0]])))*(1-(float(Single_SNP_prob[split[1]]))))*(resistance_prob)*(1-resistance_prob))
    h_dv= ((1-(float(Single_SNP_prob[split[0]])))*(1-(float(Single_SNP_prob[split[1]]))))*(1-((1-(float(Single_SNP_prob[split[0]])))*(1-(float(Single_SNP_prob[split[1]]))))*(resistance_prob)*(1-resistance_prob))
    if a_dv>0:
        a_d = math.sqrt(a_dv)
    if b_dv>0:
        b_d = math.sqrt(b_dv)
    if c_dv>0:
        c_d = math.sqrt(c_dv)
    if d_dv>0:
        d_d = math.sqrt(d_dv)
    if e_dv>0:
        e_d = math.sqrt(e_dv)
    if f_dv>0:
        f_d = math.sqrt(f_dv)
    if g_dv>0:
        g_d = math.sqrt(g_dv)
    if h_dv>0:
        h_d = math.sqrt(h_dv)
        
    if a_d >0:
        keys ["a"] = a/a_d
        keys_pair ["a"] = a/a_d
        keys_pair_single ["a"]=a/a_d
        a=a/a_d
    elif a_d==0:
        a=0
        
    if b_d >0:
        keys ["b"] = b/b_d
        keys_pair_single ["b"]=b/b_d
        b=b/b_d
    elif b_d==0:
        b=0
        
    if c_d >0:
        keys ["c"] = c/c_d
        keys_pair_single ["c"]=c/c_d
        c=c/c_d
    elif c_d==0:
        c=0
        
    if d_d >0:
        keys ["d"] = d/d_d
        d=d/d_d
    elif d_d==0:
        d=0
        
    if e_d >0:
        keys ["e"] = e/e_d
        keys_pair ["e"] = e/e_d
        keys_pair_single ["e"]=e/e_d
        e=e/e_d
    elif e_d==0:
        e=0
        
    if f_d >0:
        keys ["f"] = f/f_d
        keys_pair_single ["f"]=f/f_d
        f=f/f_d
    elif f_d==0:
        f=0
        
    if g_d >0:
        keys ["g"] = g/g_d
        keys_pair_single ["g"]=g/g_d
        g=g/g_d
    elif g_d==0:
        g=0
        
    if h_d >0:
        keys ["h"] = h/h_d
        h=h/h_d
    elif h_d==0:
        h=0

    counts = [a,b,c,d,e,f,g,h]
    counts_pair = [a,e]
    counts_pair_single =[a,b,c,e,f,g]
    r_final [element] = max(counts)
    r_pair_final [element] = max(counts_pair)
    r_pair_single_final [element] = max (counts_pair_single)
    r_finallist.append(element)
    description [element] = max(keys, key=keys.get)
    description_pair [element] = max(keys_pair, key=keys_pair.get)
    description_pair_single [element] = max(keys_pair_single, key=keys_pair_single.get)
    r_finalp.append(max(counts))
    r_pair_finalp.append(max(counts_pair))
    r_pair_single_finalp.append(max (counts_pair_single))   
    
with open (target_folder+"r_final.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in r_final.items():
        writer.writerow([key, value])
   
with open (target_folder+"r_pair_final.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in r_pair_final.items():
        writer.writerow([key, value])  
        
with open (target_folder+"r_pair_single_final.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in r_pair_single_final.items():
        writer.writerow([key, value])

with open (target_folder+"description.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in description.items():
       writer.writerow([key, value])

with open (target_folder+"description_pair.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in description_pair.items():
       writer.writerow([key, value])
       
with open (target_folder+"description_pair_single.csv", "w") as csv_file:
   writer=csv.writer(csv_file, delimiter ="\t")
   for key, value in description_pair_single.items():
       writer.writerow([key, value])
 

######################calculate r2 of SNP pairs with phenotype#################      
r2_final={}
r2_pair_final ={}
r2_pair_single_final ={}

r2_finals =[]
r2_finalp =[]
r2_pair_finals=[]
r2_pair_finalp=[]
r2_pair_single_finals=[]
r2_pair_single_finalp=[]


for element in r_finallist:
    r2_final [element] = (float(r_final[element])**2)
    r2_finals.append(element)
    r2_finalp.append(float(r_final[element])**2)
    
    r2_pair_final [element] = (float(r_pair_final[element])**2)
    r2_pair_finals.append(element)
    r2_pair_finalp.append(float(r_pair_final[element])**2)

    r2_pair_single_final [element] = (float(r_pair_single_final[element])**2)
    r2_pair_single_finals.append(element)
    r2_pair_single_finalp.append(float(r_pair_single_final[element])**2)

with open (target_folder+"r2_final1.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in r2_final.items():
        writer.writerow([key, value])
        
with open (target_folder+"r2_pair_final.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in r2_pair_final.items():
        writer.writerow([key, value])
        
with open (target_folder+"r2_pair_single_final.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in r2_pair_single_final.items():
        writer.writerow([key, value])
        
#########################Phenotype classification##############################
r2_description =[]
r2pair_description =[]
r2pair_single_description =[]
phenotype_SNP={}
r2_pair_finalpi =[]

for element in r2_finals:
    if description[element] in ["a","b", "c", "d"]:
        phenotype_SNP [element] = "Resistant"
        r2_description.append("Resistant")  
    elif description [element] in ["e","f", "g", "h"]:
        phenotype_SNP [element] = "Sensitive"
        r2_description.append("Sensitive")  
        
#for element in r2_pair_finals:
#    if description_pair[element] in ["a","b", "c", "d"]:
#        phenotype_SNP [element] = "Resistant"
#        r2pair_description.append("Resistant")  
#    elif description_pair [element] in ["e","f", "g", "h"]:
#        phenotype_SNP [element] = "Sensitive"
#        r2pair_description.append("Sensitive")  
#        
#for element in r2_pair_single_finals:
#    if description_pair_single[element] in ["a","b", "c", "d"]:
#        phenotype_SNP [element] = "Resistant"
#        r2pair_single_description.append("Resistant")  
#    elif description_pair_single [element] in ["e","f", "g", "h"]:
#        phenotype_SNP [element] = "Sensitive"
#        r2pair_single_description.append("Sensitive")  

with open (target_folder+"phenotype_SNP1.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in phenotype_SNP.items():
        writer.writerow([key, value])