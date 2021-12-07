#%%
import pandas as pd
from tqdm import tqdm
import numpy as np

def get_accessions_with_peptide_match(peptides_ifile):
    peptides = pd.read_table(peptides_ifile)
    peptides = peptides[(peptides['QValue'] <= 0.01) & (peptides['Decoy/Contaminant/Target'] =='T')]
    peptides['accs'] = peptides['Protein Accession'].str.split('|')
    accs_list = peptides['accs'].to_list()
    protein_accessions_found_in_mass_spec = []
    for accs in accs_list:
        for acc in  accs:
            protein_accessions_found_in_mass_spec.append(acc)
    protein_accessions_found_in_mass_spec = set(protein_accessions_found_in_mass_spec)
    return protein_accessions_found_in_mass_spec

def perform_subsample(orfs, subsample_size):
    subsampled_orfs = orfs.sample(weights='FL',n=subsample_size, replace=True)
    genes_found_in_ms = subsampled_orfs[subsampled_orfs['protein_found_in_ms']=='Found']['pr_gene'].unique()
    genes_not_found_in_ms = subsampled_orfs[subsampled_orfs['protein_found_in_ms']=='Not Found']['pr_gene'].unique()
    return len(genes_found_in_ms), len(genes_not_found_in_ms)
    

    
def subsample_protein_gene(peptides_ifile, orf_ifile, ofile, num_iterations=100, step_size=50000):
    orfs = pd.read_table(orf_ifile, usecols=['base_acc', 'FL', 'pr_gene'])
    protein_accessions_found_in_mass_spec = get_accessions_with_peptide_match(peptides_ifile)
    orfs['protein_found_in_ms'] = orfs['base_acc'].apply(lambda acc: 'Found' if acc in protein_accessions_found_in_mass_spec else 'Not Found')
    total = orfs['FL'].sum()

    with open(ofile, 'w') as out:
        out.write('\t'.join(['size', 'category', 'min','max','mean','sd']) + '\n')
        for subsample_size in tqdm(range(0,total,step_size)):
            array_num_genes_found = []
            array_num_genes_not_found = []
            for iteration in range(num_iterations):
                num_genes_found_in_ms, num_genes_not_found_in_ms = perform_subsample(orfs, subsample_size)
                array_num_genes_found.append(num_genes_found_in_ms)
                array_num_genes_not_found.append(num_genes_not_found_in_ms)
            array_num_genes_found = np.array(array_num_genes_found)
            array_num_genes_not_found = np.array(array_num_genes_not_found)

            out.write(f'{subsample_size}\tfound_in_ms\t{array_num_genes_found.min()}\t{array_num_genes_found.max()}\t{array_num_genes_found.mean()}\t{array_num_genes_found.std()}\n')
            out.write(f'{subsample_size}\tnot_found_in_ms\t{array_num_genes_not_found.min()}\t{array_num_genes_not_found.max()}\t{array_num_genes_not_found.mean()}\t{array_num_genes_not_found.std()}\n')
        
        
        






orf_ifile = '/Users/bj8th/Documents/Sheynkman-Lab/Data/JURKAT_06-06-2021/protein_gene_rename/jurkat_orf_refined_gene_update.tsv'
peptides_ifile = '/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics-Analysis/revisions/PacBioHybrid_NoInferenceSearch/AllPeptides.psmtsv'
ofile = './jurkat.protein_subsample.tsv'
subsample_protein_gene(peptides_ifile, orf_ifile, './test.protein_subsample.tsv',num_iterations=2, step_size=50000)


# %%
