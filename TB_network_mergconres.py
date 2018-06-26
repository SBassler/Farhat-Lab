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
from itertools import groupby, chain
import nxviz
from nxviz.plots import CircosPlot

anti ="streptomycin"
resistance_folder = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/summary_table_resistance2.tsv"
directory =  "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/"
target_folder = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/fullrun_res/"
filefolder = "/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/"
r2folder="/n/data1/hms/dbmi/farhat/sbassler/"+anti+"/files/r2/"
#antibiotic = 9 ##Rifampicin-antibiotics_dicts[18]##
#numberz=3000
########################Loading already processed strains########################
#input_file = filefolder+"filedict_"+str(numberz+500)+".csv"
#used_strains=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                used_strains.append(row[0])
##############################Antibiotics data##################################
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
gene={}
#gene=defaultdict(str)
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

input_file = filefolder+"TB_GO.csv"
genes=[]
bioprocess={}
molfunct={}
i=0
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = ";")
        for row in reader:
            if row:
                if i ==0:
                    i+=1
                elif i >=1:
                    genes.append(row[0])
                    bioprocess[row[0]] = row[2]
                    molfunct[row[0]] = row[3]
                    i+=1

#################################################################################
################for combination increase i<1001#############
#file_dict ={}
#file_list =[]
#name_list =[]
#lineages_file={}
##i=0
#for root, dir, files in os.walk(directory):
#        for name in files:
#                if name.endswith('.vcf'):
#                    name = name.split(".")
#                    names = re.findall(r"(\w+\d+)", name[0])
#                    filename = str(names[0])
#                    if (name[0]) in Counter(used_strains):
#                        if filename in antibiotics_dicts[antibiotic]:
#                            if antibiotics_dicts[antibiotic][filename] in "RS":
##                                if i < 3501:
#                                file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
#                                file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
#                                name_list.append(name[0])
#                                if (directory+name[0]+"/fast-lineage-caller/"+name[0]+".lineage"):
#                                    lineages_file [name[0]] = directory+name[0]+"/fast-lineage-caller/"+name[0]+".lineage"
##                                i +=1

#lineages={}
#lineage_list=[]
#for key, value in lineages_file.items():
#    if value:
#        try:
#            with open(value, "r") as csvfile:
#                reader = csv.reader(csvfile)
#                for row in reader:
#                    if row:
#        #                name =re.findall(r"\W(\w+)", file.split(".")[0]) [-1]
#                        if re.findall("\d\.\d", row[0]):
#                            lineages[key] = re.findall("\d\.\d", row[0])[0]
#                            lineage_list.append(re.findall("\d\.\d", row[0])[0])
#                        elif re.findall("\d", row[0]):
#                            lineages[key] = re.findall("\d", row[0])[0]
#                            lineage_list.append(re.findall("\d", row[0])[0])
#            #            lineages[name] = re.findall('\d.\d', row[1])
#            #            lineage.append(re.findall('\d.\d', row[1]))
#        except:
#            pass
##    
#with open (target_folder+"lineages_0-3500.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in lineages.items():
#        writer.writerow([key, value])
                                
#input_file = filefolder+"Single_SNP_prob_0-3500.csv"
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
#

#input_file = target_folder+"lineages_0-3500.csv"
#lineages1=[]
#lineages2=[]
#lineages3=[]
#lineages4=[]
#lineages5=[]
#lineages6=[]
#lineages7=[]
#i=0
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                try:
#                    number=re.findall("\d", row[1])[0]
#                    if number:
#                        number=re.findall("\d", row[1])[0]
#                        if i<=500:
#                            lineages1.append(int(number))
#                            i+=1
#                except:
#                    i+=1                
                
#input_file = filefolder+"list1/Multi_SNP_prob.csv"
#Multi_SNP_prob ={}
#Multi_SNP_probs =[]
#Multi_SNP_probp =[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                if float(row[1]) > 0:
#                    Multi_SNP_prob [row[0]] = float(row[1])
#                    Multi_SNP_probs.append(row[0])
#                    Multi_SNP_probp.append(float(row[1]))
                
##############################For 500 Isolates######################################
#filefolder_r = filefolder +"files/r2/"
#r2_finalf={}
#r2_finalsf =[]
#for root, dir, files in os.walk(filefolder_r):
#    for file in files:
#        with open(filefolder_r+file, "r") as csvfile:
#                reader = csv.reader(csvfile, delimiter = "\t")
#                for row in reader:
#                    if row:
#                        r2_finalf [row[0]] = float(row[1])
#                        r2_finalsf.append(row[0])
#
#        
#filefolder_d = filefolder +"files/description/"      
#phenotype_SNPf ={}
#phenos =[]
#for root, dir, files in os.walk(filefolder_d):
#    for file in files:
#        with open(filefolder_d+file, "r") as csvfile:
#                reader = csv.reader(csvfile, delimiter = "\t")
#                for row in reader:
#                    if row:
#                        phenotype_SNPf [row[0]] = row[1]
#                        phenos.append(row[0])
                                
#r2_final={}
#r2_finals =[]
#r2_finalp =[]
#phenotype_SNP ={}
#r2_description=[]
#for element in (Counter(phenos) & Counter (r2_finalsf)):
#    if float(r2_finalf[element]) >0:
#        r2_final [element] = float(r2_finalf[element])
#        r2_finals.append(element)
#        r2_finalp.append(float(r2_finalf[element]))
#        phenotype_SNP [element] = str(phenotype_SNPf [element])
#        r2_description.append(str(phenotype_SNPf [element]))
    
#with open (target_folder+"r2finals_all_1000.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in r2_final.items():
#        writer.writerow([key, value])                 
#
#with open (target_folder+"phenotype_SNP_all_1000.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in phenotype_SNP.items():
#        writer.writerow([key, value]) 
                
#######################merging description r2 resistant and sensitive##################################
#input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_500.csv"
#r2_resistant5=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row[1] == "Resistant":
#                r2_resistant5.append(row[0])
       
input_file = r2folder+"r2_final_500_0_re.csv"
r2_final5re = defaultdict(float)
r2_final5l =[]
r2_final5r = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final5l.append(row[0])
                r2_final5re[row[0]]=float(row[1])
        
#for element in (Counter(r2_resistant5)& Counter(r2_final5l)):
#                    r2_final5r [element] = r2_final5re[element]    
#                    
#del r2_resistant5
#del r2_final5re
#del r2_final5l
#                
#input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_1000.csv"
#r2_resistant10=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row[1] == "Resistant":
#                r2_resistant10.append(row[0])

input_file = r2folder+"r2_final_1000_0_re.csv"
r2_final10re = defaultdict(float)
r2_final10l =[]
r2_final10r = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final10l.append(row[0])
                r2_final10re[row[0]]=float(row[1])
        
#for element in (Counter(r2_resistant10)& Counter(r2_final10l)):
#                    r2_final10r [element] = r2_final10re[element]    
#                    
#del r2_resistant10
#del r2_final10re
#del r2_final10l
#
#input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_1500.csv"
#r2_resistant15=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row[1] == "Resistant":
#                r2_resistant15.append(row[0])

input_file = r2folder+"r2_final_1500_0_re.csv"
r2_final15re = defaultdict(float)
r2_final15l =[]
r2_final15r = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final15l.append(row[0])
                r2_final15re[row[0]]=float(row[1])
        
#for element in (Counter(r2_resistant15)& Counter(r2_final15l)):
#                    r2_final15r [element] = r2_final15re[element]    
#                    
#del r2_resistant15
#del r2_final15re
#del r2_final15l
#
#
#input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_2000.csv"
#r2_resistant20=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row[1] == "Resistant":
#                r2_resistant20.append(row[0])

input_file = r2folder+"r2_final_2000_0_re.csv"
r2_final20re = defaultdict(float)
r2_final20l =[]
r2_final20r = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final20l.append(row[0])
                r2_final20re[row[0]]=float(row[1])
        
#for element in (Counter(r2_resistant20)& Counter(r2_final20l)):
#                    r2_final20r [element] = r2_final20re[element]    
#                    
#del r2_resistant20
#del r2_final20re
#del r2_final20l
#
#input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_2500.csv"
#r2_resistant25=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row[1] == "Resistant":
#                r2_resistant25.append(row[0])

input_file = r2folder+"r2_final_2500_0_re.csv"
r2_final25re = defaultdict(float)
r2_final25l =[]
r2_final25r = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final25l.append(row[0])
                r2_final25re[row[0]]=float(row[1])
        
#for element in (Counter(r2_resistant25)& Counter(r2_final25l)):
#                    r2_final25r [element] = r2_final25re[element]    
#                    
#del r2_resistant25
#del r2_final25re
#del r2_final25l
#
#input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_3000.csv"
#r2_resistant30=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row[1] == "Resistant":
#                r2_resistant30.append(row[0])

input_file = r2folder+"r2_final_3000_0_re.csv"
r2_final30re = defaultdict(float)
r2_final30l =[]
r2_final30r = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final30l.append(row[0])
                r2_final30re[row[0]]=float(row[1])
        
#for element in (Counter(r2_resistant30)& Counter(r2_final30l)):
#                    r2_final30r [element] = r2_final30re[element]    
#                    
#del r2_resistant30
#del r2_final30re
#del r2_final30l
#
#input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_3500.csv"
#r2_resistant35=[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row[1] == "Resistant":
#                r2_resistant35.append(row[0])

input_file = r2folder+"r2_final_3500_0_re.csv"
r2_final35re = defaultdict(float)
r2_final35l =[]
r2_final35r = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                r2_final35l.append(row[0])
                r2_final35re[row[0]]=float(row[1])
        
#for element in (Counter(r2_resistant35)& Counter(r2_final35l)):
#                    r2_final35r [element] = r2_final35re[element]    
#                    
#del r2_resistant35
#del r2_final35re
#del r2_final35l

r2_final = dict(Counter(r2_final5re) + Counter(r2_final10re) + Counter(r2_final15re) + Counter(r2_final20re) + Counter(r2_final25re) + Counter(r2_final30re) + Counter(r2_final35re))

with open (target_folder+"r2_final_0-3500_re.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in dict(r2_final).items():
        writer.writerow([key, value])    
#
r2_finals=[]
r2_finalp=[]
for key, value in dict(r2_final).items():
    r2_finals.append(key)
    r2_finalp.append(value)
##
##########################Merging Single_SNP#####################################
                                
filefolder_s = filefolder +"Single_SNP_re/"
Single_SNP_prob={}
Single_SNP_counter =[]
singler=[]
first=0
second=0
for root, dir, files in os.walk(filefolder_s):
    for file in files:
        with open(filefolder_s+file, "r") as csvfile:
                reader = csv.reader(csvfile, delimiter = "\t")
                for row in reader:
                    if row:
                        if row[0] not in singler:
                            singler.append(row[0])
                            Single_SNP_prob [row[0]] = float(row[1])
                            Single_SNP_counter.append(row[0])
                        elif row[0] in singler:
                            first= Single_SNP_prob [row[0]]
                            second=(((first*(Counter(Single_SNP_counter)[row[0]]))+float(row[1]))/((Counter(Single_SNP_counter)[row[0]])+1)) 
                            Single_SNP_prob.update({row[0]: second})
                            Single_SNP_counter.append(row[0])                               
        
with open (target_folder+"Single_SNP_prob_0-3500re.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in Single_SNP_prob.items():
        writer.writerow([key, value])
      
Single_SNP_probs =[]
Single_SNP_probp =[]        
for key, value in dict(Single_SNP_prob).items():
    Single_SNP_probs.append(key)
    Single_SNP_probp.append(value)
###################################Read#########################################
        
#input_file = filefolder+"fullrun_res/"+"r2_final_0-3500_re.csv"
#r2_final ={}
#r2_finals =[]
#r2_finalp =[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                r2_final [row[0]] = float(row[1])
#                r2_finals.append(row[0])
#                r2_finalp.append(float(row[1]))
#
#input_file = filefolder+"fullrun_res/"+"Single_SNP_prob_0-3500re.csv"
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


############################Manhatten p values#################################    
#df = pd.DataFrame({"SNP":r2_finals, "pvalue":r2_finalp, "phenotype": r2_description})
#df.phenotype = df.phenotype.astype('category')
#df.phenotype = df.phenotype.cat.set_categories(['Resistant', "Sensitive"], ordered=True)
#df = df.sort_values('phenotype')
#
## How to plot gene vs. -log10(pvalue) and colour it by chromosome?
#df['r2'] = df.pvalue
#df['ind'] = range(len(df))
#df_grouped = df.groupby(('phenotype'))
#
#fig1 = plt.figure()
#ax = fig1.add_subplot(111)
#colors = ['red','green']
#x_labels = []
#x_labels_pos = []
#for num, (name, group) in enumerate(df_grouped):
#    group.plot(kind='scatter', x='ind', y='r2',color=colors[num % len(colors)], ax=ax)
#    x_labels.append(name)
#    x_labels_pos.append((group['ind'].iloc[-1   ] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
#ax.set_xticks(x_labels_pos)
#ax.set_xticklabels(x_labels)
#ax.set_xlim([0, len(df)])
#ax.set_ylim([0, 2000])
#ax.set_xlabel('phenotype')
##fig1.savefig(target_folder+'Manhatten.png')
#fig1.set_size_inches(18.5, 10.5)
#fig1.savefig(target_folder+'Manhatten_0-1000r.png', dpi=3000)
#plt.close()
####################stat learning#############################################
#clone_key ={}
#names=[]
#i=0
#resistant = 0
#sensitive = 0
#clone_lists = defaultdict()
#for key,val in file_dict.items():
#    names = re.findall(r"(\w+\d+)", key)
#    filename = str(names[0])
#    if filename in antibiotics_dicts[antibiotic]:
#        if antibiotics_dicts[antibiotic][filename] in "RS":
#            pos = []
#            with open(val) as vcffile:
#                clone_lists[i]=[record.POS for record in vcf.Reader(vcffile) if record.is_snp]
#                clone_key [i] = antibiotics_dicts[antibiotic][filename]
#                i+=1
#
##resistance_prob = (resistant/sensitive)
##with open(os.path.join(target_folder,"resistance_prob.txt"), "w") as file1:
##    toFile = str(resistance_prob)
##    file1.write(toFile)
#
#pos=[]
#snp1 =[]
#snp2=[]
#allcombination=[]
#allpos=[]
#pheno ={}
#test =0
#c=0
#scores ={}
#scorep=0
#scorec=0
#scorepl=[]
#scorecl=[]
#scoresl=[]
#label=[]
#for file in file_list:
#    if clone_key[c] == "R" or clone_key[c] == "S":
##        pos=[]
##        allcombination=[]
##        with open(file) as vcffile:
##            vcfReader = vcf.Reader(vcffile)
##            for record in vcfReader:
##                pos.append(record.POS)
##                allpos.append(record.POS)
##            i=0
##            for first in pos[:-1]:
##                snp1=first
##                for second in pos[i+1:]:
##                    snp2=second
##                    allcombination.append(str(first)+"_"+str(second))
##                i+=1
##        scorep=0
##        scorec=0
##        for element in (Counter(allcombination) & Counter (r2_finals)):
##            if (float(r2_final [element])) > 30:
##                split=element.split("_")
##                if Single_SNP_prob[split[0]] < 0.5 and Single_SNP_prob[split[1]] < 0.5:
##                    if phenotype_SNP [element] == "Resistant":
##                        scorep = scorep + (float(r2_final [element]))
##                    elif phenotype_SNP [element] == "Sensitive":
##                        scorec = scorec + (float(r2_final [element]))
##        scoresl.append(scorep+scorec)
##        scorepl.append(scorep)
##        scorecl.append(scorec)
##        scores [file] = scorep+scorec
#        if clone_key[c] == "R":
#            label.append(int(0))
#        elif clone_key[c] == "S":
#            label.append(int(1))
###       else:
###           label.append(int(2))
##        #if score >0:
##        #    pheno [file] = "Commensal"
##        #elif score <0:
##        #    pheno [file] = "Pathogenic"
##        #elif score ==0:
##        #    pheno [file] = "NA"
##        #if pheno [file] == clone_key[c]:
##        #    test +=1
#    c +=1

############################SMV
#    
## Loading some example data
##X, y = iris_data()
##X = X[:,[0, 2]]
#Xa = np.asarray(scorecl, dtype=np.float64)
#Xb = np.asarray(scorepl, dtype=np.float64)
#dataset=np.dstack([np.float64(Xa),np.float64(Xb)])
#X = dataset.reshape(len(Xa),-1)
#y=np.asarray(label)
#
#
## Initializing Classifiers
#clf1 = LogisticRegression(random_state=0)
#clf2 = RandomForestClassifier(random_state=0)
#clf3 = SVC(random_state=0, probability=True)
#eclf = EnsembleVoteClassifier(clfs=[clf1, clf2, clf3],
#                              weights=[2, 1, 1], voting='soft')
#
## Plotting Decision Regions
#gs = gridspec.GridSpec(2, 2)
#fig2 = plt.figure()#figsize=(10, 8))
#
#labels = ['Logistic Regression',
#          'Random Forest',
#          'RBF kernel SVM',
#          'Ensemble']
#
#for clf, lab, grd in zip([clf1, clf2, clf3, eclf],
#                         labels,
#                         itertools.product([0, 1],
#                         repeat=2)):
#    clf.fit(X, y)
#    ax = plt.subplot(gs[grd[0], grd[1]])
#    fig2 = plot_decision_regions(X=X, y=y,
#                                clf=clf, legend=2)
#    plt.title(lab)
#plt.xlabel("Sensitive score", size=14)
#plt.ylabel("Resistant score", size=14)
##fig2.figure.savefig(target_folder+'statistial_learning.png')
#fig2.figure.savefig(target_folder+'statistial_learning_30_0-5.pdf',format="PDF", dpi=1000)
#plt.close()
#######################Bayesion################################################
#clf1 = LogisticRegression(random_state=1)
#clf2 = RandomForestClassifier(random_state=1)
#clf3 = GaussianNB()
#clf4 = SVC()
#
#gs = gridspec.GridSpec(2, 2)
#
#fig5 = plt.figure()#figsize=(10,8))
#
#labels = ['Logistic Regression', 'Random Forest', 'Naive Bayes', 'SVM']
#for clf, lab, grd in zip([clf1, clf2, clf3, clf4],
#                         labels,
#                         itertools.product([0, 1], repeat=2)):
#
#    clf.fit(X, y)
#    ax = plt.subplot(gs[grd[0], grd[1]])
#    fig5 = plot_decision_regions(X=X, y=y, clf=clf, legend=2)
#    plt.title(lab)
#    #plt.xlabel("Commensal score", size=14)
#    #plt.ylabel("Pathogenic score", size=14)
#    
##fig5.figure.savefig(target_folder+'statistial_learning2.png')
#fig5.figure.savefig(target_folder+'statistial_learning2_30_0-5.pdf',format="PDF", dpi=1000)
#plt.close()
###############SNPs highly connected in pathogenic network#####################
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
scoredef=50
for element in r2_finals:
#    if phenotype_SNP [element] == "Resistant":
    if r2_final[element] >= scoredef:
        split=element.split("_")
        if (split[1]+"_"+split[0]) not in r2maxp:
            if Single_SNP_prob[split[0]] < 0.99 and Single_SNP_prob[split[1]] < 0.99:
                if int(split[0]) not in network_counts_p:
                    network_counts_p.append(int(split[0]))
                if int(split[1]) not in network_counts_p:
                    network_counts_p.append(int(split[1]))
                r2maxp.append(element)
                value=r2_final[element]
                if split[0] not in resistance_comb:
                    resistance_comb[split[0]] = value
                elif split[0] in resistance_comb:
                    first=resistance_comb[split[0]]
                    resistance_comb.update({split[0]: first+value})
                if split[1] not in resistance_comb:
                    resistance_comb[split[1]] = value
                elif split[1] in resistance_comb:
                    second=resistance_comb[split[0]]
                    resistance_comb.update({split[1]: second+value})
                
loc_gene=defaultdict(str)
i=0
for element in (sorted(network_counts_p)):
    position[element] = i
    for ident in names:
        if int(chrom_start [ident]) <= int(element) <= int(chrom_end[ident]):
            if gene [ident]:
                loc_gene [element] = gene [ident]
            elif not gene [ident]:
                loc_gene [element] = ident
    node_list.append(tuple((element,(loc_gene[element]),i)))
    i+=1
    
#G_p.add_nodes_from(node_list)
#r2maxpl=[]
#for element in r2_finals:
#    #if phenotype_SNP [element] == "Resistant":
#    if r2_final[element] >= 50:
#        split=element.split("_")
#        if (split[1]+"_"+split[0]) not in r2maxpl:
#            if Single_SNP_prob[split[0]] < 0.99 and Single_SNP_prob[split[1]] < 0.99:
#                G_p.add_edge((int(split[0]), loc_gene[int(split[0])], position[int(split[0])]), (int(split[1]), loc_gene[int(split[1])], position[int(split[1])]), weight=r2_final[element])
#                edges,weights = zip(*nx.get_edge_attributes(G_p,'weight').items())
#                r2maxpl.append(element)             

highest_1000={}         
highest_1000l=[]       
for element in resistance_comb:                
#        if resistance_comb[element] > 10000:
    highest_1000 [element] = resistance_comb[element]
    highest_1000l.append(resistance_comb[element])
            
    
#with open (target_folder+"highest_0-3500re_50_0-99_new.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in highest_1000.items():
#        writer.writerow([key, value]) 


#            
##results=[]
##results2=[]
##results3=[]
##results = list(map(int, network_counts_p))
##results2= sorted(results)
##results3= list(map(str, results2))
##                    
##G_p.add_nodes_from(sorted(network_counts_p))
##degrees = sorted(G_p.degree, key=lambda x: x[1], reverse=True)
##G_p.remove_nodes_from(sorted(network_counts_p))
##G_p.add_nodes_from(degrees)

#fig3 = plt.figure(figsize=(90,90))            
#nx.draw_circular(G_p,node_color='g', with_labels=True,font_size=5, edge_color=weights, edge_cmap=plt.cm.Blues, node_size=5)
#plt.axis('equal')
##plt.show()
#fig3.savefig(target_folder+'networkx_0-3500_50_0-99.pdf', format="PDF", dpi=9000)
#plt.close()
#########gene and snp score###############
node_list=[]
network_counts_i=[]
lineage_gene={}
result_list=[]
lineage_pair={}
gene_list=[]
gene_score=[]
snp_score={}
snp_list=[]
gene_score={}
#for element in r2maxp:
#    split=element.split("_")
#    if loc_gene [int(split[0])] == "embB":
#        snp_list.append(split[1])
##        if loc_gene [int(split[1])] in ["purM","purF","purH"]:
#        result_list.append(element)
#        gene_list.append(loc_gene [int(split[1])])
#        if split[1] not in snp_score:
#            snp_score [split[1]] = resistance_comb [split[1]]
#        elif split[1] in snp_score:
#            first= snp_score [split[1]]
#            snp_score.update({split[1]:first+resistance_comb [split[1]]})
##        if loc_gene [int(split[0])]:
##            if loc_gene [int(split[1])] not in gene_score:
##                gene_score [loc_gene [int(split[1])]] = resistance_comb [split[1]]
##            elif loc_gene [int(split[1])] in gene_score:
##                first=gene_score [loc_gene [int(split[1])] ]
##                gene_score.update({split[1]:first+resistance_comb [split[1]]})
#
#    elif loc_gene [int(split[1])] == "embB":
#        snp_list.append(split[0])
##        if loc_gene [int(split[0])] in ["purM","purF","purH"]:
#        result_list.append(element)
#        gene_list.append(loc_gene [int(split[0])])
#        if split[1] not in snp_score:
#            snp_score [split[0]] = resistance_comb [split[0]]
#        elif split[1] in snp_score:
#            first= snp_score [split[0]]
#            snp_score.update({split[0]:first+resistance_comb [split[0]]})
##            if loc_gene [int(split[0])]:
##                if loc_gene [int(split[0])] not in gene_score:
##                    gene_score [loc_gene [int(split[0])]] = resistance_comb [split[0]]
##                elif loc_gene [int(split[1])] in gene_score:
##                    first=gene_score [loc_gene [int(split[0])] ]
##                    gene_score.update({split[1]:first+resistance_comb [split[0]]})
##        
#        
for element in r2maxp:
    split=element.split("_")
    if split[0] not in snp_score:
        snp_score [split[0]] = resistance_comb [split[0]]
    elif split[0] in snp_score:
        first=snp_score [split[0]] 
        snp_score.update({split[0]:first+resistance_comb [split[0]]})
    if split[1] not in snp_score:
        snp_score [split[1]] = resistance_comb [split[1]]
    elif split[1] in snp_score:
        first=snp_score [split[1]] 
        snp_score.update({split[1]:first+resistance_comb [split[1]]})
    
            
for element in snp_score:
    if loc_gene [int(element)]:
        if loc_gene [int(element)] not in gene_score:
            gene_score [loc_gene [int(element)]] = snp_score[element]
        elif loc_gene [int(element)] in gene_score:
            first=gene_score [loc_gene [int(element)]] 
            gene_score.update({loc_gene [int(element)]:first+snp_score[element]})
#            
with open (target_folder+"snp_score_"+anti+"_"+str(scoredef)+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in snp_score.items():
        writer.writerow([key, value])   
        
with open (target_folder+"gene_score_"+anti+"_"+str(scoredef)+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in gene_score.items():
        writer.writerow([key, value])           
        

#bioprocess={}
bioprocess_score={}
#molfunct={}
molfunction_score={}
current_gene1=[]
#for element in network_counts_p:
for element in snp_score:
    if int(element) in loc_gene:
        if str(element) in resistance_comb:
            if loc_gene[int(element)] in bioprocess:
                current_gene1= loc_gene[int(element)]
                process1= bioprocess[current_gene1].split(",")
                process2 = [i.replace('"', '') for i in process1]
                for pro in process2:
                    if pro:
                        pro=pro.strip()
#                        if pro not in bioprocess_score:
#                            bioprocess_score[pro] =  resistance_comb[str(element)]
#                        elif pro in bioprocess_score:
#                            first= bioprocess_score[pro]
#                            bioprocess_score.update({pro: first + resistance_comb[str(element)]})
                        if pro not in bioprocess_score:
                            bioprocess_score[pro] =  snp_score[str(element)]
                        elif pro in bioprocess_score:
                            first= bioprocess_score[pro]
                            bioprocess_score.update({pro: first + snp_score[str(element)]})
#                        
current_gene2=[]
#for element in network_counts_p:
for element in snp_score:
    if int(element) in loc_gene:
        if str(element) in resistance_comb:
            if loc_gene[int(element)] in molfunct:
                current_gene2= loc_gene[int(element)]
                funct1= molfunct[current_gene2].split(",")
                funct2 = [i.replace('"', '') for i in funct1]
                for fun in funct2:
                    if fun:
                        fun=fun.strip()
#                        if fun not in molfunction_score:
#                            molfunction_score[fun] =  resistance_comb[str(element)]
#                        elif fun in molfunction_score:
#                            fun.strip()
#                            first= molfunction_score[fun]
#                            molfunction_score.update({fun: first + resistance_comb[str(element)]})
                        if fun not in molfunction_score:
                            molfunction_score[fun] =  snp_score[str(element)]
                        elif fun in molfunction_score:
                            fun.strip()
                            first= molfunction_score[fun]
                            molfunction_score.update({fun: first + snp_score[str(element)]})
        

with open (target_folder+"molfunction_score_"+anti+"_"+str(scoredef)+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in molfunction_score.items():
        writer.writerow([key, value]) 
                
with open (target_folder+"bioprocess_score_"+anti+"_"+str(scoredef)+".csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in bioprocess_score.items():
        writer.writerow([key, value]) 
##########################Circle#################################################
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
#position ={}
#node_list=[]
#for element in r2_finals:
##    if phenotype_SNP [element] == "Resistant":
#    if r2_final[element] >= 50:
#        split=element.split("_")
#        if (split[1]+"_"+split[0]) not in r2maxp:
#            if Single_SNP_prob[split[0]] < 0.9 and Single_SNP_prob[split[1]] < 0.9:
#                if int(split[0]) not in network_counts_p:
#                    network_counts_p.append(int(split[0]))
#                if int(split[1]) not in network_counts_p:
#                    network_counts_p.append(int(split[1]))
#                r2maxp.append(element)
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
#colour_dict={}
#group_dict={}
#Mb1=[]
#Mb2=[]
#Mb3=[]
#Mb4=[]
#Mb5=[]
#loc_gene=defaultdict(str)
#i=0
#gene_list=[]
##for element in range(4411532):
#for element in (sorted(network_counts_p)):
#    position[element] = i
#    for ident in names:
#        if int(chrom_start [ident]) <= int(element) <= int(chrom_end[ident]):
#            if gene [ident]:
#                loc_gene [element] = gene [ident]
#                gene_list.append(gene [ident])
#            elif not gene [ident]:
#                loc_gene [element] = ident
#                gene_list.append(ident)
##    node_list.append(tuple((element,(loc_gene[element]),i)))
#    if element <= 1000000:
#        Mb1.append(element)
#        #colour_dict [tuple((element,(loc_gene[element]),i))] = "green"
#        #group_dict [tuple((element,(loc_gene[element]),i))] = "0-1Mb"
#    elif 1000000 < element <= 2000000:
#        Mb2.append(element)
##        colour_dict [tuple((element,(loc_gene[element]),i))] = "black"
##        group_dict [tuple((element,(loc_gene[element]),i))] = "1-2Mb"
#    elif 2000000 < element <= 3000000:
#        Mb3.append(element)
##        colour_dict [tuple((element,(loc_gene[element]),i))] = "red"
##        group_dict [tuple((element,(loc_gene[element]),i))] = "2-3Mb"   
#    elif 3000000 < element <= 4000000:
#        Mb4.append(element)
##        colour_dict [tuple((element,(loc_gene[element]),i))] = "blue"
##        group_dict [tuple((element,(loc_gene[element]),i))] = "3-4Mb"
#    elif 4000000 < element:
#        Mb5.append(element)
##        colour_dict [tuple((element,(loc_gene[element]),i))] = "orange"
##        group_dict [tuple((element,(loc_gene[element]),i))] = "4-4.42"
#    i+=1
##
##bioprocess={}
#bioprocess_score={}
##molfunct={}
#molfunction_score={}
#current_gene1=[]
##for element in network_counts_p:
#for element in snp_score:
#    if int(element) in loc_gene:
#        if str(element) in resistance_comb:
#            if loc_gene[int(element)] in bioprocess:
#                current_gene1= loc_gene[int(element)]
#                process1= bioprocess[current_gene1].split(",")
#                process2 = [i.replace('"', '') for i in process1]
#                for pro in process2:
#                    if pro:
#                        pro=pro.strip()
##                        if pro not in bioprocess_score:
##                            bioprocess_score[pro] =  resistance_comb[str(element)]
##                        elif pro in bioprocess_score:
##                            first= bioprocess_score[pro]
##                            bioprocess_score.update({pro: first + resistance_comb[str(element)]})
#                        if pro not in bioprocess_score:
#                            bioprocess_score[pro] =  snp_score[str(element)]
#                        elif pro in bioprocess_score:
#                            first= bioprocess_score[pro]
#                            bioprocess_score.update({pro: first + snp_score[str(element)]})
#                        
#current_gene2=[]
##for element in network_counts_p:
#for element in snp_score:
#    if int(element) in loc_gene:
#        if str(element) in resistance_comb:
#            if loc_gene[int(element)] in molfunct:
#                current_gene2= loc_gene[int(element)]
#                funct1= molfunct[current_gene2].split(",")
#                funct2 = [i.replace('"', '') for i in funct1]
#                for fun in funct2:
#                    if fun:
#                        fun=fun.strip()
##                        if fun not in molfunction_score:
##                            molfunction_score[fun] =  resistance_comb[str(element)]
##                        elif fun in molfunction_score:
##                            fun.strip()
##                            first= molfunction_score[fun]
##                            molfunction_score.update({fun: first + resistance_comb[str(element)]})
#                        if fun not in molfunction_score:
#                            molfunction_score[fun] =  snp_score[str(element)]
#                        elif fun in molfunction_score:
#                            fun.strip()
#                            first= molfunction_score[fun]
#                            molfunction_score.update({fun: first + snp_score[str(element)]})
#        
#
#
#
#with open (target_folder+"molfunction_score_pyrazinamide_all_25.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in molfunction_score.items():
#        writer.writerow([key, value]) 
#                
#with open (target_folder+"bioprocess_score_pyrazinamide_all_25.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in bioprocess_score.items():
#        writer.writerow([key, value]) 
#
#
#G_p.add_nodes_from(Mb1, colour="green", group="0-1Mb")
#G_p.add_nodes_from(Mb2, colour="black", group="1-2Mb") 
#G_p.add_nodes_from(Mb3, colour="red", group="2-3Mb") 
#G_p.add_nodes_from(Mb4, colour="blue", group="3-4Mb") 
#G_p.add_nodes_from(Mb5, colour="orange", group="4-4.5Mb")
#r2maxpl =[]
#r2maxpv=[]
#weight_dict={}
#for element in r2_finals:
##    if phenotype_SNP [element] == "Resistant":
#    if r2_final[element] >= 50:
#        split=element.split("_")
#        if (split[1]+"_"+split[0]) not in r2maxpl:
#            if Single_SNP_prob[split[0]] < 0.99 and Single_SNP_prob[split[1]] < 0.99:
#                G_p.add_edge(int(split[0]), int(split[1]) , weight=r2_final[element])
#                #weight_dict[tuple(((int(split[0]), loc_gene[int(split[0])], position[int(split[0])]), (int(split[1]), loc_gene[int(split[1])], position[int(split[1])]))) = r2_final[element]
#                #edges,weights = zip(*nx.get_edge_attributes(G_p,'weight').items())
#                r2maxpl.append(element)
#                r2maxpv.append(r2_final[element])
#
##count=Counter(sorted(r2maxpv))
##colors = range(len(count))
##dec={}
##i=0
##for key, value in count.items():
##    dec [key] = i
##    i+=1
#
#color_node=nx.get_node_attributes(G_p,'colour')    
#for element in r2maxpl:
#    split=element.split("_")
#    #G_p [int(split[0])] [int(split[1])] ["color"] = dec[r2_final[element]]
#    G_p [int(split[0])] [int(split[1])] ["color"] = color_node [int(split[0])]
#
#fig3 = plt.figure(figsize=(30,30))
#c = CircosPlot(G_p, node_color="colour", node_grouping="group", group_order="alphabetically" ,group_label_position="middle", group_label_color=True ,edge_width="weight", edge_cmap=plt.cm.Blues, edge_color="color")
#c.draw()
##plt.show()
##fig3.savefig(target_folder+'network_nxviz_colors_single.png', format="PNG", dpi=400)
#plt.savefig(target_folder+'circos_0-3500_50_0-99.png', format="PNG", dpi=400)
#plt.close()
########################Neural SNPs##############################################
#
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
#for element in r2_finals:
#    if r2_final[element] >= 50:
#        split=element.split("_")
#        if Single_SNP_prob[split[0]] < 0.99 and Single_SNP_prob[split[1]] < 0.99:
#            if (split[1]+"_"+split[0]) not in r2maxp:
#                r2maxp.append(element)
#        
##    elif phenotype_SNP [element] == "Sensitive":
##        if r2_final[element] >= 1300:
##            split=element.split("_")
##            r2maxc.append(element)
#
#lineages_max={}
#lineage_list=[]
#matches = [[] for x in range(len(file_list))]
#r2max = r2maxp+r2maxc
#f=0
#c=0
#for key, value in file_dict.items():
#    if clone_key[c] == "R" or clone_key[c] == "S":
#        pos=[]
#        allcombination=[]
#        with open(value) as vcffile:
#            vcfReader = vcf.Reader(vcffile)
#            for record in vcfReader:
#                pos.append(record.POS)
#                allpos.append(record.POS)
#            i=0
#            for first in pos[:-1]:
#                snp1=first
#                for second in pos[i+1:]:
#                    snp2=second
#                    allcombination.append(str(first)+"_"+str(second))
#                i+=1
#        match=[]
#        for element in (Counter(allcombination) & Counter (r2maxp)):
#            if key in lineages:
#                if element not in lineages_max:
#                    lineages_max [element] = [lineages[key]]
#                elif element in lineages_max:
#                    lineage_list=lineages_max[element]
#                    lineage_list.append(lineages[key])
#                    lineages_max.update({element: lineage_list})
#            match.append(element)
#        for element in r2maxp:
#            if element in match:
#                matches[f].append(1)
#            elif element not in match:
#                matches[f].append(0)
#        f+=1
#        
##        match=[]
##        for element in (Counter(str(pos)) & Counter (network_counts_p_single)):
##            match.append(element)
##        for element in network_counts_p_single:
##            if element in match:
##                matches[f].append(1)
##            elif element not in match:
##                matches[f].append(0)    
##        f+=1
#        
#    c +=1
#
#
##with open (target_folder+"lineages_max_50_0-9.csv", "w") as csv_file:
##    writer=csv.writer(csv_file, delimiter ="\t")
##    for key, value in lineages_max.items():
##        writer.writerow([key, value])
#
#X=np.array(matches)
#y = np.array([label]).T
##
## sigmoid function
#def nonlin(x,deriv=False):
#    if(deriv==True):
#        return x*(1-x)
#    return 1/(1+np.exp(-x))
#    
##Xa = np.array(scorecl, dtype=np.float64)
##Xb = np.array(scorepl, dtype=np.float64)
##dataset=np.dstack([np.float64(Xa),np.float64(Xb)])
##X = dataset.reshape(len(Xa),-1)
#
## seed random numbers to make calculation
## deterministic (just a good practice)
#np.random.seed(1)
#
## initialize weights randomly with mean 0
#syn0 = 2*np.random.random((len(r2maxp),1)) - 1
#
#for iter in range(10000000):
#
#    # forward propagation
#    l0 = X
#    l1 = nonlin(np.dot(l0,syn0))
#
#    # how much did we miss?
#    l1_error = y - l1
#
#    # multiply how much we missed by the 
#    # slope of the sigmoid at the values in l1
#    l1_delta = l1_error * nonlin(l1,True)
#
#    # update weights
#    syn0 += np.dot(l0.T,l1_delta)  
#
#deepscores ={}
#for i in range(len(name_list)):
#    deepscores [name_list[i]] = l1[i]
#
#with open (target_folder+"deepscores_0-3500_50_0-9.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in deepscores.items():
#        writer.writerow([key, value]) 
############################2 layer neurals SNPs#################################
#deeplayers={}
#np.random.seed(1)
##
## sigmoid function
#def nonlin(x,deriv=False):
#    if(deriv==True):
#        return x*(1-x)
#    return 1/(1+np.exp(-x))
#
#i=0
## randomly initialize our weights with mean 0
#syn0 = 2*np.random.random((len(r2maxp),2)) - 1
#syn1 = 2*np.random.random((2,1)) - 1
#
#for j in range(10000000):
#
#	# Feed forward through layers 0, 1, and 2
#    l0 = X
#    l1 = nonlin(np.dot(l0,syn0))
#    l2 = nonlin(np.dot(l1,syn1))
#
#    # how much did we miss the target value?
#    l2_error = y - l2
#    
#    if (j% 1000000) == 0:
#        deeplayers [i] = np.mean(np.abs(l2_error))
#        i+=1
#        
#    # in what direction is the target value?
#    # were we really sure? if so, don't change too much.
#    l2_delta = l2_error*nonlin(l2,deriv=True)
#
#    # how much did each l1 value contribute to the l2 error (according to the weights)?
#    l1_error = l2_delta.dot(syn1.T)
#    
#    # in what direction is the target l1?
#    # were we really sure? if so, don't change too much.
#    l1_delta = l1_error * nonlin(l1,deriv=True)
#
#    syn1 += l1.T.dot(l2_delta)
#    syn0 += l0.T.dot(l1_delta)
#
#with open (target_folder+"deeplayers_0-3500_50_0-9.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in deeplayers.items():
#        writer.writerow([key, value]) 


#result_list=["1068432_840018"]
#

#location1= "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/files/r2b/"      
#value={}
#i=0
#for root, dir, files in os.walk(location1):
#        for name in files:
#            value [i] = location1+name
#            i+=1
#            
#location2= "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/files/letters/"      
#letter={}
#i=0
#for root, dir, files in os.walk(location2):
#        for name in files:
#            letter [i] = location2+name
#            i+=1
#
#
#i=0
#result_dict = {}
#snp_pairs=[]
#snp_value={}
#result_score={}
#score_list=[]
#letter_list=[]
#for i in range (70):
#    i+=1
#    snp_pairs=[]
#    snp_value={}
#    pair_letters ={}
#    with open(location1+"r2_final"+str(i)+".csv", "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                snp_pairs.append(row [0])
#                snp_value [row[0]] = float (row[1])
#    with open(location2+"description"+str(i)+".csv", "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                pair_letters [row[0]] = row[1]
#    if (Counter(result_list)& Counter(snp_pairs)):
#        for element in (Counter(result_list)& Counter(snp_pairs)):
#    #       if snp_value [element] >= 5:
#            if element not in result_dict:
#                result_dict [element] = [pair_letters[element]]
#                result_score [element] = snp_value [element]
#            elif element in result_dict:
#                letter_list=result_dict[element]
#                letter_list.append(pair_letters[element])
#                result_dict.update({element: letter_list})
#                score_list=result_score [element]
#                result_score.update({element: (score_list+(snp_value [element]))})
#
#result_dict2={}
#for key,value in result_dict.items():
#    result_dict2 [key] = Counter (value)
#    
#with open (target_folder+"SNP_scores.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in result_score.items():
#        writer.writerow([key, value]) 
#                
#with open (target_folder+"SNP_letters.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in result_dict2.items():
#        writer.writerow([key, value]) 
#
#gene_dict={}        
#for element in result_dict2:
#    split=element.split("_")
#    if (loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]) not in gene_dict:
#        gene_dict [loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]] = result_dict2[element]
#    elif (loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]) in gene_dict:
#        first= gene_dict [loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]]
#        gene_dict.update({loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]:(Counter(first)+Counter(result_dict2[element]))})
#        
#with open (target_folder+"Gene_letters.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in gene_dict.items():
#        writer.writerow([key, value]) 
            
#i=0
#pair_letters_all ={}
#letter_list=[]
#element_list=[]
#for i in range (70):
#    i+=1
##    with open(location1+"r2_final"+str(i)+".csv", "r") as csvfile:
##        reader = csv.reader(csvfile, delimiter = "\t")
##        for row in reader:
##            if row:
##                if row[0] not in snp_value:
##                    snp_pairs.append(row [0])
##                    snp_value [row[0]] = float (row[1])
##                elif row[1] in snp_value:
##                    first=snp_value [row[0]]
##                    snp_value.update({row[0]: first+float (row[1])})    
#    with open(location2+"description"+str(i)+".csv", "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                if row[0] not in pair_letters_all:
#                    element_list.append(row[0])
#                    pair_letters_all [row[0]] = [row[1]]  
#                elif row[0] in pair_letters_all:
#                    letter_list=pair_letters_all[row[0]]
#                    letter_list.append(row[1])
#                    pair_letters_all.update({row[1]: letter_list})
#   
#pair_letters_all2={}
#for key,value in pair_letters_all.items():
#    pair_letters_all2 [key] = Counter (value)
#             
#with open (target_folder+"Letters_all.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in pair_letters_all2.items():
#        writer.writerow([key, value])     

#input_file = filefolder+"letters_all.csv"
#pair_letters_all2 ={}
#element_list =[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                pair_letters_all2 [row[0]] = row[1]
#                element_list.append(row[0])
#
#SNP_letters_all={}
#for element in (Counter(result_list)& Counter(element_list)):
#    SNP_letters_all [element] = pair_letters_all2 [element]
#    
#
#with open (target_folder+"SNP_letters_all_new_50.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in SNP_letters_all.items():
#        writer.writerow([key, value]) 
#                
#gene_dict_all={}        
#for element in SNP_letters_all:
#    split=element.split("_")
#    if (loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]) not in gene_dict_all:
#        gene_dict_all [loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]] = pair_letters_all2[element]
#    elif (loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]) in gene_dict_all:
#        first= gene_dict_all [loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]]
#        gene_dict_all.update({loc_gene [int(split[0])]+"_"+loc_gene [int(split[1])]:first+"_"+pair_letters_all2[element]})
#        
#with open (target_folder+"Gene_letters_all_new_50.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in gene_dict_all.items():
#        writer.writerow([key, value])  

#                
#r2_final_500 = defaultdict()
#r2_final_1000 = defaultdict()
#r2_final_1500 = defaultdict()
#r2_final_2000 = defaultdict()
#r2_final_2500 = defaultdict()
#r2_final_3000 = defaultdict()
#r2_final_3500 = defaultdict()                
##                
#number=30
#for count in range(number):
#    input_file = r2folder+"r2_final_500_"+str(count)+"_re.csv"
#    r2_final5re = defaultdict(float)
#    r2_final5l =[]
#    r2_final5r = defaultdict(float)
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    r2_final5l.append(row[0])
#                    r2_final5re[row[0]]=float(row[1])
#            
#    #for element in (Counter(r2_resistant5)& Counter(r2_final5l)):
#    #                    r2_final5r [element] = r2_final5re[element]    
#    #                    
#    #del r2_resistant5
#    #del r2_final5re
#    #del r2_final5l
#    #                
#    #input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_1000.csv"
#    #r2_resistant10=[]
#    #with open(input_file, "r") as csvfile:
#    #        reader = csv.reader(csvfile, delimiter = "\t")
#    #        for row in reader:
#    #            if row[1] == "Resistant":
#    #                r2_resistant10.append(row[0])
#    
#    input_file = r2folder+"r2_final_1000_"+str(count)+"_re.csv"
#    r2_final10re = defaultdict(float)
#    r2_final10l =[]
#    r2_final10r = defaultdict(float)
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    r2_final10l.append(row[0])
#                    r2_final10re[row[0]]=float(row[1])
#            
#    #for element in (Counter(r2_resistant10)& Counter(r2_final10l)):
#    #                    r2_final10r [element] = r2_final10re[element]    
#    #                    
#    #del r2_resistant10
#    #del r2_final10re
#    #del r2_final10l
#    
#    #input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_1500.csv"
#    #r2_resistant15=[]
#    #with open(input_file, "r") as csvfile:
#    #        reader = csv.reader(csvfile, delimiter = "\t")
#    #        for row in reader:
#    #            if row[1] == "Resistant":
#    #                r2_resistant15.append(row[0])
#    
#    input_file = r2folder+"r2_final_1500_"+str(count)+"_re.csv"
#    r2_final15re = defaultdict(float)
#    r2_final15l =[]
#    r2_final15r = defaultdict(float)
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    r2_final15l.append(row[0])
#                    r2_final15re[row[0]]=float(row[1])
#            
#    #for element in (Counter(r2_resistant15)& Counter(r2_final15l)):
#    #                    r2_final15r [element] = r2_final15re[element]    
#    #                    
#    #del r2_resistant15
#    #del r2_final15re
#    #del r2_final15l
#    #
#    #
#    #input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_2000.csv"
#    #r2_resistant20=[]
#    #with open(input_file, "r") as csvfile:
#    #        reader = csv.reader(csvfile, delimiter = "\t")
#    #        for row in reader:
#    #            if row[1] == "Resistant":
#    #                r2_resistant20.append(row[0])
#    #
#    input_file = r2folder+"r2_final_2000_"+str(count)+"_re.csv"
#    r2_final20re = defaultdict(float)
#    r2_final20l =[]
#    r2_final20r = defaultdict(float)
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    r2_final20l.append(row[0])
#                    r2_final20re[row[0]]=float(row[1])
#    #        
#    #for element in (Counter(r2_resistant20)& Counter(r2_final20l)):
#    #                    r2_final20r [element] = r2_final20re[element]    
#    #                    
#    #del r2_resistant20
#    #del r2_final20re
#    #del r2_final20l
#    #
#    #input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_2500.csv"
#    #r2_resistant25=[]
#    #with open(input_file, "r") as csvfile:
#    #        reader = csv.reader(csvfile, delimiter = "\t")
#    #        for row in reader:
#    #            if row[1] == "Resistant":
#    #                r2_resistant25.append(row[0])
#    
#    input_file = r2folder+"r2_final_2500_"+str(count)+"_re.csv"
#    r2_final25re = defaultdict(float)
#    r2_final25l =[]
#    r2_final25r = defaultdict(float)
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    r2_final25l.append(row[0])
#                    r2_final25re[row[0]]=float(row[1])
#    #        
#    #for element in (Counter(r2_resistant25)& Counter(r2_final25l)):
#    #                    r2_final25r [element] = r2_final25re[element]    
#    #                    
#    #del r2_resistant25
#    #del r2_final25re
#    #del r2_final25l
#    #
#    #input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_3000.csv"
#    #r2_resistant30=[]
#    #with open(input_file, "r") as csvfile:
#    #        reader = csv.reader(csvfile, delimiter = "\t")
#    #        for row in reader:
#    #            if row[1] == "Resistant":
#    #                r2_resistant30.append(row[0])
#    #
#    input_file = r2folder+"r2_final_3000_"+str(count)+"_re.csv"
#    r2_final30re = defaultdict(float)
#    r2_final30l =[]
#    r2_final30r = defaultdict(float)
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    r2_final30l.append(row[0])
#                    r2_final30re[row[0]]=float(row[1])
#    #        
#    #for element in (Counter(r2_resistant30)& Counter(r2_final30l)):
#    #                    r2_final30r [element] = r2_final30re[element]    
#    #                    
#    #del r2_resistant30
#    #del r2_final30re
#    #del r2_final30l
#    #
#    #input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_3500.csv"
#    #r2_resistant35=[]
#    #with open(input_file, "r") as csvfile:
#    #        reader = csv.reader(csvfile, delimiter = "\t")
#    #        for row in reader:
#    #            if row[1] == "Resistant":
#    #                r2_resistant35.append(row[0])
#    #
#    input_file = r2folder+"r2_final_3500_"+str(count)+"_re.csv"
#    r2_final35re = defaultdict(float)
#    r2_final35l =[]
#    r2_final35r = defaultdict(float)
#    with open(input_file, "r") as csvfile:
#            reader = csv.reader(csvfile, delimiter = "\t")
#            for row in reader:
#                if row:
#                    r2_final35l.append(row[0])
#                    r2_final35re[row[0]]=float(row[1])
#    #        
#    #for element in (Counter(r2_resistant35)& Counter(r2_final35l)):
#    #                    r2_final35r [element] = r2_final35re[element]    
#    #                    
#    #del r2_resistant35
#    #del r2_final35re
#    #del r2_final35l
#    #
#    #r2_final = dict(Counter(r2_final5r) + Counter(r2_final10r) + Counter(r2_final15r))# + Counter(r2_final20r) + Counter(r2_final25r) + Counter(r2_final30r) + Counter(r2_final35r))
#    r2_final = dict(Counter(r2_final5re) + Counter(r2_final10re) + Counter(r2_final15re) + Counter(r2_final20re) + Counter(r2_final25re) + Counter(r2_final30re) + Counter(r2_final35re))
#    r2_final_500 [count]= Counter(r2_final5re)
#    r2_final_1000 [count]= Counter(r2_final10re)
#    r2_final_1500 [count]= Counter(r2_final15re)
#    r2_final_2000 [count]= Counter(r2_final20re)
#    r2_final_2500 [count]= Counter(r2_final25re)
#    r2_final_3000 [count]= Counter(r2_final30re)
#    r2_final_3500 [count]= Counter(r2_final35re)
##    with open (filefolder+"collection/"+"r2_final_0-3500r_"+str(count)+".csv", "w") as csv_file:
##        writer=csv.writer(csv_file, delimiter ="\t")
##        for key2, value2 in dict(r2_final).items():
##            writer.writerow([key2, value2])    
#
#
#results=[]
#pdist=defaultdict()
#p09999=defaultdict()
#p0999=defaultdict()
#p099=defaultdict()
#max_count=defaultdict()
#pelement=[]
#i=0
#for a in range (number):
##    for b in range (number):
##        for c in range (number):
##            for d in range (number):
##                for e in range (number):
##                    for f in range (number):
##                        for g in range (number):
#    merge= r2_final_500 [a] + r2_final_1000 [a] + r2_final_1500 [a] + r2_final_2000 [a] + r2_final_2500 [a] + r2_final_3000 [a] + r2_final_3500 [a]
#    pelement=[]
#    for key, value in dict(merge).items():
#       pelement.append(float(value))
##                            pdist [i] = pelement
#    test=np.array(pelement)
#    test2=test.round(decimals=2)
#    p09999[i]= np.percentile(test2, 99.99)
#    p0999[i]= np.percentile(test2, 99.9)
#    p099[i]= np.percentile(test2, 99)
#    max_count[i] = (max(test2))
#    i+=1
#
#with open (target_folder+"p09999_"+str(number)+".csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in p09999.items():
#        writer.writerow([key, value]) 
#
#with open (target_folder+"p0999_"+str(number)+".csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in p0999.items():
#        writer.writerow([key, value]) 
#        
#with open (target_folder+"p099_"+str(number)+".csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in p099.items():
#        writer.writerow([key, value]) 
#
#with open (target_folder+"max_count_"+str(number)+".csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in max_count.items():
#        writer.writerow([key, value]) 
