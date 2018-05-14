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
target_folder = "/n/data1/hms/dbmi/farhat/sbassler/rifampicin/files/fullrun_counter/"
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

################for combination increase i<1001#############
file_dict ={}
file_list =[]
name_list =[]
lineages_file=[]
i=0
for root, dir, files in os.walk(directory):
        for name in files:
                if name.endswith('.vcf'):
                    name = name.split(".")
                    names = re.findall(r"(\w+\d+)", name[0])
                    filename = str(names[0])
                    if filename in antibiotics_dicts[antibiotic]:
                        if antibiotics_dicts[antibiotic][filename] in "RS":
                            if i < 1001:
                                file_list.append(directory+name[0]+"/pilon/"+name[0]+".vcf")
                                file_dict [name[0]] = directory+name[0]+"/pilon/"+name[0]+".vcf"
                                name_list.append(name[0])
                                lineages_file.append(directory+name[0]+"/fast-lineage-caller/"+name[0]+".lineage")
                                i +=1

#lineages={}
#lineage=[]
#for file in lineages_file:
#    try:
#        with open(file, "r") as csvfile:
#            reader = csv.reader(csvfile)
#            for row in reader:
#                name =re.findall(r"\W(\w+)", file.split(".")[0]) [-1]
#                lineages[file] = re.findall("\d", row[0])
#                lineage.append(re.findall('\d', row[0]))
#    #            lineages[name] = re.findall('\d.\d', row[1])
#    #            lineage.append(re.findall('\d.\d', row[1]))
#    except:
#        pass
#    
#with open (target_folder+"lineages_0-1000.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in lineages.items():
#        writer.writerow([key, value])
                                
input_file = filefolder+"Single_SNP_prob_0-1000.csv"
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
#
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
#filefolder_r = filefolder +"r2/"
#r2_finalf={}
#r2_finalsf =[]
#for root, dir, files in os.walk(filefolder_r):
#    for file in files[0:10]:
#        with open(filefolder_r+file, "r") as csvfile:
#                reader = csv.reader(csvfile, delimiter = "\t")
#                for row in reader:
#                    if row:
#                        r2_finalf [row[0]] = float(row[1])
#                        r2_finalsf.append(row[0])
#        
#filefolder_d = filefolder +"description/"      
#phenotype_SNPf ={}
#phenos =[]
#for root, dir, files in os.walk(filefolder_d):
#    for file in files[0:10]:
#        with open(filefolder_d+file, "r") as csvfile:
#                reader = csv.reader(csvfile, delimiter = "\t")
#                for row in reader:
#                    if row:
#                        phenotype_SNPf [row[0]] = row[1]
#                        phenos.append(row[0])
#                                
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
#    
#with open (filefolder +"r2m/"+"r2finals_all_1500.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in r2_final.items():
#        writer.writerow([key, value])                 
#
#with open (filefolder +"description2m/"+"phenotype_SNP_all_1500.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in phenotype_SNP.items():
#        writer.writerow([key, value]) 

########################Using default dict#####################################

input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_500.csv"
r2_resistant5=[]
r2_sensitive5=[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row[1] == "Resistant":
                r2_resistant5.append(row[0])
            elif row[1] == "Sensitive":
                r2_sensitive5.append(row[0])
                
input_file = filefolder +"descriptionm/"+"phenotype_SNP_all_1000.csv"
r2_resistant10=[]
r2_sensitive10=[]
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row[1] == "Resistant":
                r2_resistant10.append(row[0])
            elif row[1] == "Sensitive":
                r2_sensitive10.append(row[0])

input_file = filefolder +"r2m/"+"r2finals_all_500.csv"
r2_final5r = defaultdict(float)
r2_final5s = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                if row[0] in Counter (r2_resistant5):
                    r2_final5r [row[0]] = float(row[1])
                elif row[0] in Counter (r2_sensitive5):
                    r2_final5s [row[0]] = float(row[1])
                
input_file = filefolder +"r2m/"+"r2finals_all_1000.csv"
r2_final10r = defaultdict(float)
r2_final10s = defaultdict(float)
with open(input_file, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")
        for row in reader:
            if row:
                if row[0] in Counter (r2_resistant5):
                    r2_final10r [row[0]] = float(row[1])
                elif row[0] in Counter (r2_sensitive5):
                    r2_final10s [row[0]] = float(row[1])

#r2_descriptionf=defaultdict(str)
phenotype_SNPf ={}
r2_finalr = Counter(r2_final5r) + Counter(r2_final10r)
r2_finals = Counter(r2_final5s) + Counter(r2_final10s)
for element in r2_finalr:
    phenotype_SNPf [element] = "Resistant"
for element in r2_finals:
    if element in phenotype_SNPf:
        pass
    else:
        phenotype_SNPf [element] = "Sensitive"

r2_final= dict(Counter(r2_finalr) | Counter (r2_finals))
for element in Counter(r2_finalr) & Counter(r2_finals):
    if r2_finalr [element] > r2_finals [element]:
        phenotype_SNPf.update({element: "Resistant"})
    elif r2_finals [element] > r2_finalr [element]:
        phenotype_SNPf.update({element: "Sensitive"})    

r2_finals=[]
r2_finalp=[]
phenotype_SNP={}
r2_description=[]        
for key, value in dict(r2_final).items():
    r2_finals.append(key)
    r2_finalp.append(value)
    phenotype_SNP [key] = phenotype_SNPf [key]
    r2_description.append(phenotype_SNPf [key])

with open (target_folder+"r2_final_0-1000.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in dict(r2_final).items():
        writer.writerow([key, value]) 
        
with open (target_folder+"phenotype_SNP_0-1000.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in dict(phenotype_SNP).items():
        writer.writerow([key, value])

#######################merging description r2 resistant and sensitive##################################
#filefolder_d = filefolder +"descriptionm/"     
#phenos =[]
#single=[]
#phenotype_SNPf={}
#for root, dir, files in os.walk(filefolder_d):
#    for file in files:
#        with open(filefolder_d+file, "r") as csvfile:
#                reader = csv.reader(csvfile, delimiter = "\t")
#                for row in reader:
#                    if row:
#                        if row[0] not in single:
#                            single.append(row[0])
#                            phenotype_SNPf [row[0]] = row[1]
#                            phenos.append(row[1])
#                        elif row[0] in single:
#                            phenos.append(row[1])
#
#filefolder_r = filefolder +"r2m/"
#r2_finalf={}
#r2_finalsf =[]
#singler=[]
#first=0
#second=0
#i=0
#for root, dir, files in os.walk(filefolder_r):
#    for file in files:
#        with open(filefolder_r+file, "r") as csvfile:
#                reader = csv.reader(csvfile, delimiter = "\t")
#                for row in reader:
#                    if row:
#                        if row[0] not in singler:
#                            singler.append(row[0])
#                            r2_finalf [row[0]] = float(row[1])
#                            r2_finalsf.append(row[0])
#                        elif row[0] in singler:
#                            if phenotype_SNPf [row[0]] == phenos[i]:
#                                first= r2_finalf [row[0]]
#                                second=first+float(row[1])
#                                #second=(((first*(Counter(r2_finalsf)[row[0]]))+float(row[1]))/((Counter(r2_finalsf)[row[0]])+1)) 
#                                r2_finalf.update({row[0]: second})
#                                r2_finalsf.append(row[0])
#                            else:
#                                if float(row[1])>first:
#                                   phenotype_SNPf.update({row[0]: phenos[i]})
#                                   r2_finalf.update({row[0]: float(row[1])})
#                                   r2_finalsf.append(row[0])
#                    i+=1
#                
#r2_final={}
#r2_finals =[]
#r2_finalp =[]
#phenotype_SNP ={}
#r2_description=[]
#for element in (Counter(single) & Counter (singler)):
#    if float(r2_finalf[element]) >0:
#        r2_final [element] = float(r2_finalf[element])
#        r2_finals.append(element)
#        r2_finalp.append(float(r2_finalf[element]))
#        phenotype_SNP [element] = str(phenotype_SNPf [element])
#        r2_description.append(str(phenotype_SNPf [element]))
#    
#with open (target_folder+"r2finals_all_0-1000.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in r2_final.items():
#        writer.writerow([key, value])                 
#
#with open (target_folder+"phenotype_SNP_all_0-1000.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in phenotype_SNP.items():
#        writer.writerow([key, value]) 
#        
        
########################Merging Single_SNP#####################################
#input_file = filefolder+"list1/Single_SNP_prob.csv"
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
                                
#filefolder_s = filefolder +"files/Single_SNP/"
#Single_SNP_prob={}
#Single_SNP_counter =[]
#singler=[]
#first=0
#second=0
#for root, dir, files in os.walk(filefolder_r):
#    for file in files:
#        with open(filefolder_r+file, "r") as csvfile:
#                reader = csv.reader(csvfile, delimiter = "\t")
#                for row in reader:
#                    if row:
#                        if row[0] not in singler:
#                            singler.append(row[0])
#                            Single_SNP_prob [row[0]] = float(row[1])
#                            Single_SNP_counter.append(row[0])
#                        elif row[0] in singler:
#                            first= r2_final [row[0]]
#                            second=(((first*(Counter(Single_SNP_counter)[row[0]]))+float(row[1]))/((Counter(Single_SNP_counter)[row[0]])+1)) 
#                            Single_SNP_prob.update({row[0]: second})
#                            Single_SNP_counter.append(row[0])                               
#        
#with open (target_folder+"Single_SNP_prob_0-1000.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in Single_SNP_prob.items():
#        writer.writerow([key, value])              

######Read description and r2###################################################
        
#input_file = target_folder+"r2finals_0-1000.csv"
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
#input_file = target_folder+"phenotype_SNP_0-1000.csv"
#phenotype_SNP ={}
#r2_description =[]
#with open(input_file, "r") as csvfile:
#        reader = csv.reader(csvfile, delimiter = "\t")
#        for row in reader:
#            if row:
#                phenotype_SNP [row[0]] = row[1]
#                r2_description.append(row[1])


############################Manhatten p values#################################    
df = pd.DataFrame({"SNP":r2_finals, "pvalue":r2_finalp, "phenotype": r2_description})
df.phenotype = df.phenotype.astype('category')
df.phenotype = df.phenotype.cat.set_categories(['Resistant', "Sensitive"], ordered=True)
df = df.sort_values('phenotype')

# How to plot gene vs. -log10(pvalue) and colour it by chromosome?
df['r2'] = df.pvalue
df['ind'] = range(len(df))
df_grouped = df.groupby(('phenotype'))

fig1 = plt.figure()
ax = fig1.add_subplot(111)
colors = ['red','green']
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='r2',color=colors[num % len(colors)], ax=ax)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1   ] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
ax.set_ylim([0, 2000])
ax.set_xlabel('phenotype')
#fig1.savefig(target_folder+'Manhatten.png')
fig1.set_size_inches(20, 20)
fig1.savefig(target_folder+'Manhatten_0-1000.pdf', format="PDF", dpi=5000)
plt.close()
#####################stat learning#############################################

clone_key ={}
names=[]
i=0
resistant = 0
sensitive = 0
clone_lists = defaultdict()
for key,val in file_dict.items():
    names = re.findall(r"(\w+\d+)", key)
    filename = str(names[0])
    if filename in antibiotics_dicts[antibiotic]:
        if antibiotics_dicts[antibiotic][filename] in "RS":
            if antibiotics_dicts[antibiotic][filename] == "R":
                resistant +=1
            elif antibiotics_dicts[antibiotic][filename] == "S":
                sensitive +=1
            pos = []
            with open(val) as vcffile:
                clone_lists[i]=[record.POS for record in vcf.Reader(vcffile) if record.is_snp]
                clone_key [i] = antibiotics_dicts[antibiotic][filename]
                i+=1

resistance_prob = (resistant/sensitive)
with open(os.path.join(target_folder,"resistance_prob_0-1000.txt"), "w") as file1:
    toFile = str(resistance_prob)
    file1.write(toFile)

pos=[]
snp1 =[]
snp2=[]
allcombination=[]
allpos=[]
pheno ={}
test =0
c=0
scores ={}
scorep=0
scorec=0
scorepl=[]
scorecl=[]
scoresl=[]
label=[]
for file in file_list:
    if clone_key[c] == "R" or clone_key[c] == "S":
        pos=[]
        allcombination=[]
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
        scorep=0
        scorec=0
        for element in (Counter(allcombination) & Counter (r2_finals)):
            if (float(r2_final [element])) > 30:
                split=element.split("_")
                if Single_SNP_prob[split[0]] < 0.7 and Single_SNP_prob[split[1]] < 0.7:
                    if phenotype_SNP [element] == "Resistant":
                        scorep = scorep + (float(r2_final [element]))
                    elif phenotype_SNP [element] == "Sensitive":
                        scorec = scorec + (float(r2_final [element]))
#        scoresl.append(scorep+scorec)
#        scorepl.append(scorep)
#        scorecl.append(scorec)
#        scores [file] = scorep+scorec
        if clone_key[c] == "R":
            label.append(int(0))
        elif clone_key[c] == "S":
            label.append(int(1))
##       else:
##           label.append(int(2))
#        #if score >0:
#        #    pheno [file] = "Commensal"
#        #elif score <0:
#        #    pheno [file] = "Pathogenic"
#        #elif score ==0:
#        #    pheno [file] = "NA"
#        #if pheno [file] == clone_key[c]:
#        #    test +=1
    c +=1

############################SMV
#    
## Loading some example data
##X, y = iris_data()
##X = X[:,[0, 2]]
Xa = np.asarray(scorecl, dtype=np.float64)
Xb = np.asarray(scorepl, dtype=np.float64)
dataset=np.dstack([np.float64(Xa),np.float64(Xb)])
X = dataset.reshape(len(Xa),-1)
y=np.asarray(label)


# Initializing Classifiers
clf1 = LogisticRegression(random_state=0)
clf2 = RandomForestClassifier(random_state=0)
clf3 = SVC(random_state=0, probability=True)
eclf = EnsembleVoteClassifier(clfs=[clf1, clf2, clf3],
                              weights=[2, 1, 1], voting='soft')

# Plotting Decision Regions
gs = gridspec.GridSpec(2, 2)
fig2 = plt.figure()#figsize=(10, 8))

labels = ['Logistic Regression',
          'Random Forest',
          'RBF kernel SVM',
          'Ensemble']

for clf, lab, grd in zip([clf1, clf2, clf3, eclf],
                         labels,
                         itertools.product([0, 1],
                         repeat=2)):
    clf.fit(X, y)
    ax = plt.subplot(gs[grd[0], grd[1]])
    fig2 = plot_decision_regions(X=X, y=y,
                                clf=clf, legend=2)
    plt.title(lab)
plt.xlabel("Sensitive score", size=14)
plt.ylabel("Resistant score", size=14)
#fig2.figure.savefig(target_folder+'statistial_learning.png')
fig2.figure.savefig(target_folder+'statistial_learning_0-1000_30_0-7.pdf',format="PDF", dpi=3000)
plt.close()
######################Bayesion################################################
clf1 = LogisticRegression(random_state=1)
clf2 = RandomForestClassifier(random_state=1)
clf3 = GaussianNB()
clf4 = SVC()

gs = gridspec.GridSpec(2, 2)

fig5 = plt.figure()#figsize=(10,8))

labels = ['Logistic Regression', 'Random Forest', 'Naive Bayes', 'SVM']
for clf, lab, grd in zip([clf1, clf2, clf3, clf4],
                         labels,
                         itertools.product([0, 1], repeat=2)):

    clf.fit(X, y)
    ax = plt.subplot(gs[grd[0], grd[1]])
    fig5 = plot_decision_regions(X=X, y=y, clf=clf, legend=2)
    plt.title(lab)
    #plt.xlabel("Commensal score", size=14)
    #plt.ylabel("Pathogenic score", size=14)
    
#fig5.figure.savefig(target_folder+'statistial_learning2.png')
fig5.figure.savefig(target_folder+'statistial_learning2_0-1000_30_0-7.pdf',format="PDF", dpi=3000)
plt.close()
################SNPs highly connected in pathogenic network#####################
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
G_p.add_nodes_from(node_list)
r2maxpl=[]
for element in r2_finals:
    if phenotype_SNP [element] == "Resistant":
        if r2_final[element] >= 50:
            split=element.split("_")
            if (split[1]+"_"+split[0]) not in r2maxpl:
                if Single_SNP_prob[split[0]] < 0.9 and Single_SNP_prob[split[1]] < 0.9:
                    G_p.add_edge((int(split[0]), loc_gene[int(split[0])], position[int(split[0])]), (int(split[1]), loc_gene[int(split[1])], position[int(split[1])]), weight=r2_final[element])
                    edges,weights = zip(*nx.get_edge_attributes(G_p,'weight').items())
                    r2maxpl.append(element)
                    

#highest_1000={}         
#highest_1000l=[]       
#for element in resistance_comb:                
#        if resistance_comb[element] > 10000:
#            highest_1000 [element] = resistance_comb[element]
#            highest_1000l.append(resistance_comb[element])
            
    
#with open (target_folder+"highest_1000.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in highest_1000.items():
#        writer.writerow([key, value]) 
            
#results=[]
#results2=[]
#results3=[]
#results = list(map(int, network_counts_p))
#results2= sorted(results)
#results3= list(map(str, results2))
#                    
#G_p.add_nodes_from(sorted(network_counts_p))
#degrees = sorted(G_p.degree, key=lambda x: x[1], reverse=True)
#G_p.remove_nodes_from(sorted(network_counts_p))
#G_p.add_nodes_from(degrees)

fig3 = plt.figure(figsize=(90,90))            
nx.draw_circular(G_p,node_color='g', with_labels=True,font_size=5, edge_color=weights, edge_cmap=plt.cm.Blues, node_size=5)
plt.axis('equal')
#plt.show()
fig3.savefig(target_folder+'network_p_finalsorted_0-9_0-1000_50_genes_regions.pdf', format="PDF", dpi=9000)
plt.close()
#########################Neural SNPs##############################################

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
for element in r2_finals:
    if phenotype_SNP [element] == "Resistant":
        if r2_final[element] >= 30:
            split=element.split("_")
            if Single_SNP_prob[split[0]] < 0.7 and Single_SNP_prob[split[1]] < 0.7:
                if (split[1]+"_"+split[0]) not in r2maxp:
                    r2maxp.append(element)
        
#    elif phenotype_SNP [element] == "Sensitive":
#        if r2_final[element] >= 1300:
#            split=element.split("_")
#            r2maxc.append(element)

matches = [[] for x in range(len(file_list))]
r2max = r2maxp+r2maxc
f=0
c=0
for file in file_list:
    if clone_key[c] == "R" or clone_key[c] == "S":
        pos=[]
        allcombination=[]
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
        match=[]
        for element in (Counter(allcombination) & Counter (r2maxp)):
            match.append(element)
        for element in r2maxp:
            if element in match:
                matches[f].append(1)
            elif element not in match:
                matches[f].append(0)
        f+=1
        
#        match=[]
#        for element in (Counter(str(pos)) & Counter (network_counts_p_single)):
#            match.append(element)
#        for element in network_counts_p_single:
#            if element in match:
#                matches[f].append(1)
#            elif element not in match:
#                matches[f].append(0)    
#        f+=1
        
    c +=1
#
X=np.array(matches)
y = np.array([label]).T

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
#for iter in range(1000000):
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
#with open (target_folder+"deepscores.csv", "w") as csv_file:
#    writer=csv.writer(csv_file, delimiter ="\t")
#    for key, value in deepscores.items():
#        writer.writerow([key, value]) 
###########################2 layer neurals SNPs#################################
deeplayers={}
np.random.seed(1)
#
# sigmoid function
def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))

i=0
# randomly initialize our weights with mean 0
syn0 = 2*np.random.random((len(r2maxp),2)) - 1
syn1 = 2*np.random.random((2,1)) - 1

for j in range(1000000):

	# Feed forward through layers 0, 1, and 2
    l0 = X
    l1 = nonlin(np.dot(l0,syn0))
    l2 = nonlin(np.dot(l1,syn1))

    # how much did we miss the target value?
    l2_error = y - l2
    
    if (j% 100000) == 0:
        deeplayers [i] = np.mean(np.abs(l2_error))
        i+=1
        
    # in what direction is the target value?
    # were we really sure? if so, don't change too much.
    l2_delta = l2_error*nonlin(l2,deriv=True)

    # how much did each l1 value contribute to the l2 error (according to the weights)?
    l1_error = l2_delta.dot(syn1.T)
    
    # in what direction is the target l1?
    # were we really sure? if so, don't change too much.
    l1_delta = l1_error * nonlin(l1,deriv=True)

    syn1 += l1.T.dot(l2_delta)
    syn0 += l0.T.dot(l1_delta)

with open (target_folder+"deeplayers_0-1000_30_0-7.csv", "w") as csv_file:
    writer=csv.writer(csv_file, delimiter ="\t")
    for key, value in deeplayers.items():
        writer.writerow([key, value]) 