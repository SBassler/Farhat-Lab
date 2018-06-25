# Farhat-Lab
network_matrix_1-7 calculate SNP pair scores for every SNP pair correlated with resistance to a given antibiotic (example streptomycin).
Each of the scripts is calculating the SNP pair values for 500 isolates and uses the list of already used strains from the previous script to just take unused strains.
Mergconres takes all the scores, adds them and calculates SNP scores, gene scores and performs gene ontology searches.
It also has tools to search for SNPs and genes related to a specific gene or to plot SNP pairs as a network with edges representing scores of a pair or regions of the first SNP. Additionally, it can call lineages, in which a SNP pair is found. It can also perform a permutation analysis or look for letters representing a SNP pair in each subset (a-h).
