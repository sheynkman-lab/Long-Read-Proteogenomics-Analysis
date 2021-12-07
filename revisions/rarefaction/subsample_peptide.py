#%%
from os import access
import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import random
def map_peptides_to_orfs(peptides):
    orf_peptides_map = defaultdict(set)
    for _, row in peptides.iterrows():
        for acc in row['accs']:
            orf_peptides_map[acc].add(row['Full Sequence'])
    return orf_peptides_map


def perform_rarefaction(rarefaction_list, orf_peptides_map, subsample_size):
    accessions_with_peptides = set(orf_peptides_map.keys())
    subsampled_accessions = set(random.sample(rarefaction_list, subsample_size))
    subsampled_accessions_with_peptides = subsampled_accessions.intersection(accessions_with_peptides)
    peptides_subsampled = [orf_peptides_map.get(acc) for acc in subsampled_accessions_with_peptides]
    peptides_subsampled = set().union(*peptides_subsampled)
    return len(peptides_subsampled)
    
def read_peptides_file(peptides_ifile):
    peptides = pd.read_table(peptides_ifile)
    peptides = peptides[(peptides['QValue'] <= 0.01) & (peptides['Decoy/Contaminant/Target'] =='T')]
    peptides['accs'] = peptides['Protein Accession'].str.split('|')
    return peptides

def generate_rarefaction_list(orfs):
    rarefaction_list = []
    for _, row in orfs.iterrows():
        acc_list = [row['base_acc']]*row['FL']
        rarefaction_list.extend(acc_list)
    return rarefaction_list

    
def generate_peptide_rarefaction_table(peptides_ifile, orf_ifile, ofile, num_iterations=100, step_size=50000):
    orfs = pd.read_table(orf_ifile, usecols=['base_acc', 'FL', 'pr_gene'])
    orf_peptides_map = map_peptides_to_orfs(read_peptides_file(peptides_ifile))
    total = orfs['FL'].sum()
    rarefaction_list = generate_rarefaction_list(orfs)
    with open(ofile, 'w') as out:
        out.write('\t'.join(['size', 'min','max','mean','sd']) + '\n')
        peptides_found_per_subsample = defaultdict(list)
        for subsample_size in tqdm(range(0,total,step_size)):
            for iteration in range(num_iterations):
                num_peptides_found = perform_rarefaction(rarefaction_list, orf_peptides_map, subsample_size)
                peptides_found_per_subsample[subsample_size].append(num_peptides_found)
            peptides_found_per_subsample[subsample_size] = np.array(peptides_found_per_subsample[subsample_size])
            num_peptides_found_in_subsample_size = peptides_found_per_subsample[subsample_size]
            out.write(f'{subsample_size}\t{num_peptides_found_in_subsample_size.min()}\t{num_peptides_found_in_subsample_size.max()}\t{num_peptides_found_in_subsample_size.mean()}\t{num_peptides_found_in_subsample_size.std()}\n')
       
        
#%%
orf_ifile = '/Users/bj8th/Documents/Sheynkman-Lab/Data/JURKAT_06-06-2021/protein_gene_rename/jurkat_orf_refined_gene_update.tsv'
peptides_ifile = './AllPeptides_filteredDB.psmtsv'
num_iterations = 20
step_size = 25000

generate_peptide_rarefaction_table(peptides_ifile, orf_ifile, f'jurkat.filtered.peptides_subsampled.step_{step_size}.iter_{num_iterations}.tsv', num_iterations=num_iterations, step_size=step_size)





# %%
