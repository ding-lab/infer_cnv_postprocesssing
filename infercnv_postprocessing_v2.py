#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys

# input example:
# infercnv_postprocessing.py
# argv[1] = infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt
# argv[2] = infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.references.txt

# make a gene to band dictionary
karyo_map = {};
karyo_file = open("Karyotype_genes_hg38.txt");
karyo_file.readline();
for line in karyo_file:
    gene,chrom,band_spec,start,stop = line.strip().split("\t");
    if band_spec:
        band_collapsed = chrom+band_spec[0];
    else:
        band_collapsed = chrom;
    karyo_map[gene] = band_collapsed;
karyo_file.close();

# process tumor/observations output
cell_cnv_dict = {};
cell_cnv_dict_gene_level = {};
cnv_file = open(sys.argv[1]);
cells = cnv_file.readline().strip().replace('"', '').split(" ");
for cell in cells:
    cell_cnv_dict[cell] = {};
    cell_cnv_dict_gene_level[cell] = {};
for line in cnv_file:
    line = line.strip().split(" ");
    cnv_values = [float(j) for j in line[1:]];
    gene = line[0].replace('"', '');
    for i in range(len(cells)):
        cell_cnv_dict_gene_level[cells[i]][gene] = cnv_values[i];
    if gene not in karyo_map:
        continue;
    band = karyo_map[gene];
    if "p" not in band and "q" not in band:
        continue;
    for i in range(len(cells)):
        if band not in cell_cnv_dict[cells[i]]:
            cell_cnv_dict[cells[i]][band] = [];
        cell_cnv_dict[cells[i]][band].append(cnv_values[i]);
cnv_file.close();

# process non-tumor/reference output
cnv_file = open(sys.argv[2]);
cells = cnv_file.readline().strip().replace('"', '').split(" ");
for cell in cells:
    cell_cnv_dict[cell] = {};
    cell_cnv_dict_gene_level[cell] = {};
for line in cnv_file:
    line = line.strip().split(" ");
    cnv_values = [float(j) for j in line[1:]];
    gene = line[0].replace('"', '');
    for i in range(len(cells)):
        cell_cnv_dict_gene_level[cells[i]][gene] = cnv_values[i];
    if gene not in karyo_map:
        continue;
    band = karyo_map[gene];
    if "p" not in band and "q" not in band:
        continue;
    for i in range(len(cells)):
        if band not in cell_cnv_dict[cells[i]]:
            cell_cnv_dict[cells[i]][band] = [];
        cell_cnv_dict[cells[i]][band].append(cnv_values[i]);
cnv_file.close();
                
# make annotations at the arm level
for cell in cell_cnv_dict:
    for band in cell_cnv_dict[cell]:
        cnv = np.mean(cell_cnv_dict[cell][band]);
        if cnv <= 0.25:
            cnv = "0x";
        elif cnv <= 0.75:
            cnv = "0.5x";
        elif cnv <= 1.25:
            cnv = "1x";
        elif cnv <= 1.75:
            cnv = "1.5x";
        elif cnv <= 2.25:
            cnv = "2x";
        else:
            cnv = "3x";
        cell_cnv_dict[cell][band] = cnv;

# make annotations at the gene level
for cell in cell_cnv_dict_gene_level:
    for gene in cell_cnv_dict_gene_level[cell]:
        cnv = cell_cnv_dict_gene_level[cell][gene];
        if cnv <= 0.25:
            cnv = "0x";
        elif cnv <= 0.75:
            cnv = "0.5x";
        elif cnv <= 1.25:
            cnv = "1x";
        elif cnv <= 1.75:
            cnv = "1.5x";
        elif cnv <= 2.25:
            cnv = "2x";
        else:
            cnv = "3x";
        cell_cnv_dict_gene_level[cell][gene] = cnv;

# write out the tables
output = pd.DataFrame().from_dict(cell_cnv_dict,orient="index").reset_index();
output.to_csv(sys.argv[3]+"_arm_level.tsv",sep = "\t",index=False);
output2 = pd.DataFrame().from_dict(cell_cnv_dict_gene_level,orient="index").reset_index();
output2.to_csv(sys.argv[3]+"_gene_level.tsv",sep = "\t",index=False);
