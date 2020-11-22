#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from Bio import SeqIO

__author__ = 'Matthew Onorato'
__version__ = '0.0.1'
__date__ = 'May 7, 2020'
__maintainer__ = 'Matthew Onorato'
__email__ = 'monorato2@gmail.com' > 'monorato4398@sdsu.edu'
__status__ = 'final'


"""
DESCRIPTION: This script extracts locus tags from GBK files, finds the union of all 'Rv' genes, and then assigns a '1'
             if the gene is present in a given genome or a '0' if the gene is absent.
             --> This analysis does not include novel genes (e.g. those prefixed with 'MTB' or some other ID) as I'm not
                 sure how consistent the labeling was. If the labeling was inconsistent between genomes, this can lead
                 to mis-groupings for features, which would mess with the classification predictions.
             --> note: All genomes were sequenced on PacBio's RSII platform.
             
IMPORTANT: This script is not meant to be run as it uses files with private information. I wanted to include it for
           completeness.  
"""

gbks = [f for f in os.listdir('genomes/') if not f.startswith('.')]  # ensures no hidden files get included

strain_to_gene = {}

for gbk in gbks:
    with open(f'genomes/{gbk}', 'r') as handle:
        try:
            record = SeqIO.read(handle, 'genbank')
        except ValueError:
            for record in SeqIO.parse(handle, 'genbank'):
                print(record.id)

        locus_tags = []
        for feature in record.features:
            if feature.type == 'CDS':
                for qualifier, value in feature.qualifiers.items():
                    if qualifier == 'gene':
                        gene = value[0]
                    elif qualifier == 'locus_tag':
                        locus_tag = value[0]

                        # finds all 'Rv' labeled genes, but does not capture novel genes
                        if 'Rv' in locus_tag:
                            locus_tags.append(locus_tag)
                        elif 'Rv' in gene:
                            locus_tags.append(gene)
                        # else:
                        #     locus_tags.append(gene)  # for novel genes (some annotated w/ gene name, e.g. catD)

    strain_name = gbk.split('.')[0]
    strain_to_gene[strain_name] = locus_tags

# finds union of all genes, but also filters out duplicate genes (OK for this analysis)
unique_genes = set()
for vs in strain_to_gene.values():
    for v in vs:
        unique_genes.add(v)

strain_names = [gbk.split('.')[0] for gbk in gbks]
df_strains = pd.DataFrame(columns=strain_names, index=unique_genes)

# assigns '1' if gene present in strain, otherwise assigns '0'
for strain_name, strain_genes in strain_to_gene.items():
    for strain_gene in strain_genes:
        df_strains.loc[strain_gene, strain_name] = 1

df_strains.replace(np.nan, 0, inplace=True)
df_strains = df_strains.T

# df_strains.to_csv('Mtb_pangenome_analysis.csv', mode='w')


# pangenome genes from lab --> Mtb_pangenome_analysis.csv still needs strain_name replacement and DR_status
df_lab = pd.read_csv('Mtb_pangenome_analysis.csv')
df_lab.rename(columns={'Unnamed: 0': 'Isolate'}, inplace=True)

# maps generic strain ID to hide actual strain names
name_map = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in open('strain_map.txt', 'r')}
df_lab['Isolate'] = df_lab['Isolate'].map(name_map)
df_lab.set_index('Isolate', inplace=True)


df_pan = pd.read_excel('Mtb_pangenome_analysis.xlsx', usecols='A,U')
df_pan.set_index('Isolate', inplace=True)

# assigns 0 for mono-resistant strains, 1 for MDR, and 2 for XDR (note: only 4 drug susceptible strains, so not used)
df_pan['Resistance Classification'] = df_pan['Resistance Classification'].apply(
    lambda x: '0' if x == 'mono' else ('1' if x == 'MDR' or x == 'preXDR' else ('2' if x == 'XDR' else 'NA')))

strain_to_status = {idx: row['Resistance Classification'] for idx, row in df_pan.iterrows()}

# appends DR_status column (our label) to the gene presence/absence dataframe
for strain_name in df_lab.index.values:
    try:
        df_lab.loc[strain_name, 'DR_status'] = strain_to_status[strain_name]
    except KeyError:
        continue
        # df_lab.loc[strain_name, 'DR_status'] = 'panS'  # mainly for H37Rv & CDC1551 (PATRIC database says panS)

# df_lab.to_csv('Mtb_pangenome_analysis.csv', mode='w')
