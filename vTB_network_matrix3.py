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
from collections import defaultdict
from collections import OrderedDict
import random

output = "1500"
anti="streptomycin"
antibiotic = 19 ##Rifampicin-antibiotics_dicts[18], pyraziamide [16]## [6] ethambutol### Cipro [3]
numberz=(int(output)-500)

resistance_folder = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/summary_table_resistance2.tsv"
directory =  "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/"
target_folder = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/list1_"+output+"/"
filefolder = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/list1_"+output+"/"
r2folder="/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/r2/"


######################Loading already processed strains########################
#input_file = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/filedict_3500.csv"
#used_strains=[]
#used_dict={}
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                used_strains.append(row[0])
#                used_dict [row[0]] = row [1]


input_file = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/list1_"+str(numberz)+"/filedict_"+str(numberz)+".csv"
used_strains=[]
used_dict={}
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                used_strains.append(row[0])
                used_dict [row[0]] = row [1]

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
test={}
i=0
for root, dir, files in os.walk(directory):
    for name in files:
        if name.endswith('.vcf'):
            name = name.split(".")
            names = re.findall(r"(\w+\d+)", name[0])
            filename = str(names[0])
#            if (name[0]) in Counter(used_strains[numberz:(numberz+500)]):
            if (name[0]) not in Counter(used_strains):
                if filename in antibiotics_dicts[antibiotic]:
                    if antibiotics_dicts[antibiotic][filename] == "R":
                        if i <=500:
                            file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
                            file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
                            name_list.append(name)
                            i +=1
                    elif antibiotics_dicts[antibiotic][filename] == "S":
                        if i <=500:
                            file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
                            file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
                            name_list.append(name)
                            i +=1
                 
with open (target_folder+"filedict_"+output+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in used_dict.items():
        writer.writerow([key, value])
    for key, value in file_dict.items():
        writer.writerow([key, value])               
        
                    
########Single_SNP_probability and Multi_SNP_probability (1 hour)############
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

Single_SNP_prob ={}
Single_SNP_proba ={}
Single_SNP_probs =[]
Single_SNP_probp =[]
counts = Counter(allpos)
Single_SNP_proba=dict(counts)
for element in Single_SNP_proba:
    if (float(Single_SNP_proba[element])/len(file_list)) >=0.01:
        Single_SNP_prob[str(element)] = (float(Single_SNP_proba[element])/len(file_list))
        Single_SNP_probs.append(str(element))
        Single_SNP_probp.append(float(Single_SNP_proba[element])/len(file_list))
    
with open (target_folder+"Single_SNP_prob_"+output+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in Single_SNP_prob.items():
        writer.writerow([key, value])
        
with open ("/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/"+"Single_SNP_re/"+"Single_SNP_prob_"+output+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in Single_SNP_prob.items():
        writer.writerow([key, value])

for file in file_list:
    pos=[]
    with open(file) as vcffile:
        vcfReader = vcf.Reader(vcffile)
        for record in vcfReader:
            pos.append(record.POS)
            allpos.append(record.POS)
        i=0
        for first in pos[:-1]:
            if str(first) in Single_SNP_prob:
                snp1=first
                for second in pos[i+1:]:
                    if str(second) in Single_SNP_prob:
                        snp2=second
                        allcombination.append(str(first)+"_"+str(second))
                i+=1

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
    
with open (target_folder+"Multi_SNP_prob_"+output+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in Multi_SNP_prob.items():
        writer.writerow([key, value])
        
##########################Load SNPs############################################
#input_file = filefolder+"Single_SNP_prob_"+output+".csv"
#Single_SNP_prob ={}
#Single_SNP_probs =[]
#Single_SNP_probp =[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                Single_SNP_prob [row[0]] = float(row[1])
#                Single_SNP_probs.append(row[0])
#                Single_SNP_probp.append(float(row[1]))        
                
                
input_file = filefolder+"Multi_SNP_prob_"+output+".csv"
Multi_SNP_prob ={}
Multi_SNP_probs =[]
Multi_SNP_probp =[]
element_dict={}
element_list=[]
single=[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                Multi_SNP_prob [row[0]] = float(row[1])
                Multi_SNP_probs.append(row[0])
                Multi_SNP_probp.append(float(row[1]))
                split=row[0].split("_")
                if split[0] not in element_dict:
                    element_dict [split[0]] = [row[0]]
                elif split[0] in element_dict:
#                    if row[0] not in Counter(single):
#                        single.append(row[0])
#                        single.append(split[1]+"_"+split[0])
                    element_list= element_dict [split[0]]
                    element_list.append(row[0])
                    element_dict.update({split[0]: element_list})
                if split[1] not in element_dict:
                    element_dict [split[1]] = [row[0]]
                elif split[1] in element_dict:
#                    if row[0] not in Counter(single):
#                        single.append(row[0])
#                        single.append(split[1]+"_"+split[0])
                    element_list= element_dict [split[1]]
                    element_list.append(row[0])
                    element_dict.update({split[1]: element_list})
                
##################Phenotype correlations for different drugs###################
res_data={}
clone_key_first =[]
names=[]
i=0
resistant = 0
sensitive = 0
clone_names=[]
clone_lists = defaultdict()
for key,val in file_dict.items():
    names = re.findall(r"(\w+\d+)", key)
    filename = str(names[0])
    if filename in antibiotics_dicts[antibiotic]:
        if antibiotics_dicts[antibiotic][filename] in "RS":
            if antibiotics_dicts[antibiotic][filename] == "R":
                resistant +=1
                res=1
            elif antibiotics_dicts[antibiotic][filename] == "S":
                sensitive +=1
                res=0
            pos = []
            with open(val) as vcffile:
                clone_lists[key]=[record.POS for record in vcf.Reader(vcffile) if record.is_snp]
                clone_key_first.append(res)
                clone_names.append(key)
                res_data[key] = res
                i+=1

resistance_prob = (resistant/sensitive)
clone_key=np.array(clone_key_first)
#with open (target_folder+"res_data.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in res_data.items():
#        writer.writerow([key, value])

Single_SNP_prob_list=defaultdict()      
SNP_list=[] 
for element in Single_SNP_probs:
    SNP_list=[]
    for file in clone_lists:
        if int(element) in Counter(clone_lists[file]):
            SNP_list.append(1)
        elif int(element) not in Counter(clone_lists[file]):
            SNP_list.append(0)
    Single_SNP_prob_list [element] = np.array(SNP_list)
    
#with open (target_folder+"Single_SNP_prob_arrays.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in Single_SNP_prob_list.items():
#        writer.writerow([key, value])


for counter in range(1):
#    random.shuffle(clone_key)
    count_list_a =defaultdict()
    count_list_b= defaultdict()
    count_list_c= defaultdict()
    count_list_d= defaultdict()
    count_list_e= defaultdict()
    count_list_f= defaultdict()
    count_list_g = defaultdict()
    count_list_h= defaultdict()
    count=0
    length=(len(Multi_SNP_probs))
        
    #############################Through Multi_SNP_prob############################
    array_list=[]
    names=[]
    result_list=[]
    final_array = defaultdict()
    result_list=[]
    for key, value in element_dict.items():
        snp2=[]
        names=[]
        snp1 = Single_SNP_prob_list[key]
        for element in value:
            split=element.split("_")
            if split[0] == key:
                snp2.append(Single_SNP_prob_list[split[1]])
                names.append(split[1])
            elif split[1] == key:
                snp2.append(Single_SNP_prob_list[split[0]])
                names.append(split[0])
        final_array [key] = (snp2+(snp1<<1)+(clone_key<<2))
        for m in range (len(final_array[key])):
            values = Counter (final_array[key][m])
            count_list_a [key+"_"+names[m]] = values[7] /len(clone_key)
            count_list_b [key+"_"+names[m]] = values[6] /len(clone_key)
            count_list_c [key+"_"+names[m]] = values[5] /len(clone_key)
            count_list_d [key+"_"+names[m]] = values[4] /len(clone_key)
            count_list_e [key+"_"+names[m]] = values[3] /len(clone_key)
            count_list_f [key+"_"+names[m]] = values[2] /len(clone_key)
            count_list_g [key+"_"+names[m]] = values[1] /len(clone_key)
            count_list_h [key+"_"+names[m]] = values[0] /len(clone_key)
        del final_array [key]
        
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
    a_d=0
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
    for element in Multi_SNP_probs:
        ar=0
        br=0
        cr=0
        dr=0
        er=0
        fr=0
        gr=0
        hr=0
        split=element.split("_")
        counts = []
    #    counts_pair = []
    #    counts_pair_single = []
        keys ={}
    #    keys_pair={}
    #    keys_pair_single={}
        a = count_list_a[element]-(Multi_SNP_prob[element])*(resistance_prob)
        b = count_list_b[element]-((float(Single_SNP_prob[split[0]])*(1-(float(Single_SNP_prob[split[1]]))))*(resistance_prob))
        c = count_list_c[element]-((float(Single_SNP_prob[split[1]])*(1-(float(Single_SNP_prob[split[0]]))))*(resistance_prob))
        d = count_list_d[element]-(1-(Multi_SNP_prob[element]))*(resistance_prob)
        e = count_list_e[element]-(Multi_SNP_prob[element])*(1-resistance_prob)
        f = count_list_f[element]-((float(Single_SNP_prob[split[0]])*(1-(float(Single_SNP_prob[split[1]]))))*(1-resistance_prob))
        g = count_list_g[element]-((float(Single_SNP_prob[split[1]])*(1-(float(Single_SNP_prob[split[0]]))))*(1-resistance_prob))
        h = count_list_h[element]-(1-(Multi_SNP_prob[element]))*(1-resistance_prob)
        a_dv=(float(Multi_SNP_prob[element]))*(1-(float(Multi_SNP_prob[element])))*(resistance_prob)*(1-resistance_prob)
        b_dv=((float(Single_SNP_prob[split[0]]))*(1-(float(Single_SNP_prob[split[1]]))))*(1-((float(Single_SNP_prob[split[0]]))*(1-(float(Single_SNP_prob[split[1]])))))*(resistance_prob)*(1-resistance_prob)
        c_dv=((float(Single_SNP_prob[split[1]]))*(1-(float(Single_SNP_prob[split[0]]))))*(1-((float(Single_SNP_prob[split[1]]))*(1-(float(Single_SNP_prob[split[0]])))))*(resistance_prob)*(1-resistance_prob)
        d_dv=(1-(float(Multi_SNP_prob[element])))*(1-(1-(float(Multi_SNP_prob[element]))))*(resistance_prob)*(1-resistance_prob)
        e_dv=(float(Multi_SNP_prob[element]))*(1-(float(Multi_SNP_prob[element])))*(resistance_prob)*(1-resistance_prob)
        f_dv=((float(Single_SNP_prob[split[0]]))*(1-(float(Single_SNP_prob[split[1]]))))*(1-((float(Single_SNP_prob[split[0]]))*(1-(float(Single_SNP_prob[split[1]])))))*(resistance_prob)*(1-resistance_prob)
        g_dv=((float(Single_SNP_prob[split[1]]))*(1-(float(Single_SNP_prob[split[0]]))))*(1-((float(Single_SNP_prob[split[1]]))*(1-(float(Single_SNP_prob[split[0]])))))*(resistance_prob)*(1-resistance_prob)
        h_dv=(1-(float(Multi_SNP_prob[element])))*(1-(1-(float(Multi_SNP_prob[element]))))*(resistance_prob)*(1-resistance_prob)
        
        if a_dv>0:
            a_d = math.sqrt(a_dv)
            keys ["a"] = a/a_d
    #        keys_pair ["a"] = a/a_d
    #        keys_pair_single ["a"]=a/a_d
            ar=a/a_d
            
        if b_dv>0:
            b_d = math.sqrt(b_dv)
            keys ["b"] = b/b_d
    #        keys_pair_single ["b"]=b/b_d
            br=b/b_d
    
        if c_dv>0:
            c_d = math.sqrt(c_dv)
            keys ["c"] = c/c_d
    #        keys_pair_single ["c"]=c/c_d
            cr=c/c_d
            
        if d_dv>0:
            d_d = math.sqrt(d_dv)
            keys ["d"] = d/d_d
            dr=d/d_d
            
        if e_dv>0:
            e_d = math.sqrt(e_dv)
            keys ["e"] = e/e_d
    #        keys_pair ["e"] = e/e_d
    #        keys_pair_single ["e"]=e/e_d
            er=e/e_d
            
        if f_dv>0:
            f_d = math.sqrt(f_dv)        
            keys ["f"] = f/f_d
    #        keys_pair_single ["f"]=f/f_d
            fr=f/f_d
            
        if g_dv >0:
            g_d = math.sqrt(g_dv)
            keys ["g"] = g/g_d
    #        keys_pair_single ["g"]=g/g_d
            gr=g/g_d
            
        if h_dv>0:
            h_d = math.sqrt(h_dv)        
            keys ["h"] = h/h_d
            hr=h/h_d
    
        counts = [ar,br,cr,dr,er,fr,gr,hr]
    #    counts = [ar,br,cr,dr]
    #    counts_pair = [ar,er]
    #    counts_pair_single =[ar,br,cr,er,fr,gr]
        if keys:
            if max(keys, key=keys.get) in "abcd":
                r_final [element] = max(counts)
            #    r_pair_final [element] = max(counts_pair)
            #    r_pair_single_final [element] = max (counts_pair_single)
                r_finallist.append(element)
                description [element] = max(keys, key=keys.get)
            #    description_pair [element] = max(keys_pair, key=keys_pair.get)
            #    description_pair_single [element] = max(keys_pair_single, key=keys_pair_single.get)
                r_finalp.append(max(counts))
            #    r_pair_finalp.append(max(counts_pair))
            #    r_pair_single_finalp.append(max (counts_pair_single))
        
    
    
    with open (target_folder+"description_"+output+"_re.csv", "w") as csv_file:
       writer=csv.writer(csv_file, delimiter ="\t")
       for key2, value2 in description.items():
           writer.writerow([key2, value2])
     
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
    
    with open (r2folder+"r2_final_"+output+"_"+str(counter)+"_re.csv", "w") as csv_file:
        writer=csv.writer(csv_file, delimiter ="\t")
        for key3, value3 in r2_final.items():
            writer.writerow([key3, value3])