# CSS844_mod3
Code for module 3 of CSS 844.

Coexpression analysis of Wisconsin diversity panel.

filter_data_wgcna.r - takes gene expression data (TPMs) and filters according to varience, removing genes and samples that have too few entries. Data output is RData file which can be used in downstream WGCNA analysis.

module_sft.r - takes filtered gene expression data in the form of a RData file (produced by filter_data_wgcna.r) and determines the best soft thresholding power to use to create a coexpression network. Produces two graphs showing scale free model fits for each soft thresholding power tested, and mean connectivity for each tested power. 

module_detection.r - takes filtered gene expression data in the form of a RData file (produced by filter_data_wgcna.r) and soft thresholding power determined by module_sft.r. Outputs a network in the form of an adjacency matrix, cacluated in part using the selected soft thresholding power.

get_module_genes.r - creates a file listing each module and the number of genes in that module as well as individual .csv files for each module containing all the gene names in each module. 

module_traits.r - takes module data produced by module_detection.r and external phenotype data and correlates module eigengene values to external trait data. Outputs a correlation matrix that can be input into heatmap.py for visualization.

heatmap.py - takes module/trait correlation matrix and outputs a heatmap visualizing correlations between modules and traits.
