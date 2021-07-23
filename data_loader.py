#%%
from shutil import SameFileError
import pandas as pd
from Bio import SeqIO
import config
import re

sample_name = 'jurkat'
gene_isoname_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/reference_tables/gene_isoname.tsv'
pb_gene_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/protein_classification/{sample_name}_genes.tsv'
uniprot_fasta_file = f'{config.REFERENCE_DIRECTORY}/uniprot_reviewed_canonical_and_isoform.fasta'
# PSM file locations
gencode_psm_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/gencode/search_results/Task1SearchTask/AllPSMs.psmtsv'
uniprot_psm_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/uniprot/search_results/Task1SearchTask/AllPSMs.psmtsv'
pacbio_filtered_psm_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/filtered/search_results/Task1SearchTask/AllPSMs.psmtsv'
pacbio_hybrid_psm_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/hybrid/search_results/Task1SearchTask/AllPSMs.psmtsv'
# Peptide file locations
gencode_peptide_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/gencode/search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv'
uniprot_peptide_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/uniprot/search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv'
pacbio_filtered_peptide_file = f'{config.PIPELINE_RESULTS_DIRECTORY}//metamorpheus/pacbio/filtered/search_results/Task1SearchTask/AllPeptides.{sample_name}.filtered.psmtsv'
pacbio_hybrid_peptide_file = f'{config.PIPELINE_RESULTS_DIRECTORY}//metamorpheus/pacbio/hybrid/search_results/Task1SearchTask/AllPeptides.{sample_name}.hybrid.psmtsv'
# Protein group file locations
gencode_protein_group_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/gencode/search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv'
uniprot_protein_group_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/uniprot/search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv'
pacbio_filtered_protein_group_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/filtered/search_results/Task1SearchTask/AllQuantifiedProteinGroups.{sample_name}.filtered.tsv'
pacbio_hybrid_protein_group_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/hybrid/search_results/Task1SearchTask/AllQuantifiedProteinGroups.{sample_name}.hybrid.tsv'
pacbio_rescue_resolve_protein_group_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/metamorpheus/pacbio/rescue_resolve/search_results/Task1SearchTask/AllQuantifiedProteinGroups.{sample_name}.rescue_resolve.tsv'
high_confidence_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/hybrid_protein_database/{sample_name}_refined_high_confidence.tsv'
high_confidence_genes_file = f'{config.PIPELINE_RESULTS_DIRECTORY}/hybrid_protein_database/jurkat_high_confidence_genes.tsv'


# list of genes in the high confidence space, for filtering of input data
high_confidence = pd.read_table(high_confidence_file)
accs_in_hiconf_space = set(high_confidence['base_acc'].to_list())
# genes_in_hiconf_space = set(high_confidence['gene'].to_list())
genes_in_hiconf_space = set(pd.read_table(high_confidence_genes_file)['GENE'].to_list())

def read_metamorpheus_file(metamorpheus_file, gene_map):
    metamorpheus_table = pd.read_table(metamorpheus_file)
    metamorpheus_table = metamorpheus_table[metamorpheus_table['Decoy/Contaminant/Target']=='T']
    metamorpheus_table = metamorpheus_table[metamorpheus_table['QValue']<=0.01]
    metamorpheus_table['accs'] = metamorpheus_table['Protein Accession'].str.split('|')
    metamorpheus_table['genes'] = metamorpheus_table['accs'].apply(lambda accs: [gene_map[x] for x in accs])
    return metamorpheus_table

def read_peptide_file(peptides_file, gene_map):
    return read_metamorpheus_file(peptides_file, gene_map)

def read_psm_file(psm_file, gene_map):
    return read_metamorpheus_file(psm_file, gene_map)

def get_genes_in_protein_group(protein_group, gene_map):
    genes = []
    for acc in protein_group.split('|'):
        if acc in gene_map.keys():
            genes.append(gene_map[acc])
    return genes

def is_high_confidence(list_to_check, high_confidence_list):
    for x in list_to_check:
        if x in high_confidence_list:
            return True
    return False

def read_protein_group_file(protein_groups_filename, gene_map):
    protein_groups = (
        pd.read_table(
            protein_groups_filename,
            index_col=False,
            usecols=[0,1,5,6,7,8,9,16,19], 
            names=['protein_group','pg_gene','num_proteins_in_group','unique_peptides','shared_peptides','number_of_peptides','number_of_unique_peptides','dct','qval'], 
            skiprows=1)
        .query('qval <= 0.01 and dct == "T"')
        .assign(has_unique_peptides=lambda x: ~pd.isna(x['unique_peptides']))
        .assign(number_of_shared_peptides=lambda x: x['number_of_peptides'] - x['number_of_unique_peptides'])
    )
    protein_groups['accs'] = protein_groups['protein_group'].str.split('|')
    protein_groups['genes']= protein_groups['protein_group'].apply(lambda pgroup: get_genes_in_protein_group(pgroup, gene_map))
    
    return protein_groups


#***************************************************
# gene maps
#***************************************************
gencode_gene_map = pd.read_table(gene_isoname_file, names=['gene','isoname'])
gencode_gene_map = pd.Series(gencode_gene_map.gene.values, index=gencode_gene_map.isoname).to_dict()
pacbio_gene_map = pd.read_table(pb_gene_file)
pacbio_gene_map = pd.Series(pacbio_gene_map.pr_gene.values, index=pacbio_gene_map.pb).to_dict()
hybrid_gene_map = {**pacbio_gene_map, **gencode_gene_map}

uniprot_gene_map = {}
for record in SeqIO.parse(uniprot_fasta_file, format='fasta'):
    isoform = record.id.split('|')[1]
    gene_regex = re.search(r'GN=(.*)', record.description)
    if gene_regex is not None:
        gene = gene_regex.group(1).split(' ')[0]
        uniprot_gene_map[isoform] = gene
    else:
        uniprot_gene_map[isoform] = None

#***************************************************
# Metamorpheus PSMs
#***************************************************
class PSM:
    def __init__(self) -> None:
        self.gencode_psm = read_psm_file(gencode_psm_file, gencode_gene_map)
        self.gencode_psm['is_high_confidence'] = self.gencode_psm['genes'].apply(lambda genes: is_high_confidence(genes, genes_in_hiconf_space))

        self.uniprot_psm = read_psm_file(uniprot_psm_file, uniprot_gene_map)
        self.uniprot_psm['is_high_confidence'] = self.uniprot_psm['genes'].apply(lambda genes: is_high_confidence(genes, genes_in_hiconf_space))

        # pacbio_filtered_psm = read_psm_file(pacbio_filtered_psm_file)
        self.pacbio_hybrid_psm = read_psm_file(pacbio_hybrid_psm_file, hybrid_gene_map)
        self.pacbio_hybrid_psm['is_high_confidence'] = self.pacbio_hybrid_psm['accs'].apply(lambda accs: is_high_confidence(accs, accs_in_hiconf_space))

#***************************************************
# Metamorpheus Peptides
#***************************************************
class Peptide:
    def __init__(self) -> None:
        self.gencode_peptide = read_peptide_file(gencode_peptide_file, gencode_gene_map)
        self.gencode_peptide['is_high_confidence'] = self.gencode_peptide['genes'].apply(lambda genes: is_high_confidence(genes, genes_in_hiconf_space))
        # self.gencode_peptide_set = set(self.gencode_peptide['Base Sequence'].unique())

        self.uniprot_peptide = read_peptide_file(uniprot_peptide_file, uniprot_gene_map)
        self.uniprot_peptide['is_high_confidence'] = self.uniprot_peptide['genes'].apply(lambda genes: is_high_confidence(genes, genes_in_hiconf_space))
        # self.uniprot_peptide_set = set(self.uniprot_peptide['Base Sequence'].unique())
        # pacbio_filtered_peptide = read_peptide_file(pacbio_filtered_peptide_file)
        self.pacbio_hybrid_peptide = read_peptide_file(pacbio_hybrid_peptide_file, hybrid_gene_map)
        self.pacbio_hybrid_peptide['is_high_confidence'] = self.pacbio_hybrid_peptide['accs'].apply(lambda accs: is_high_confidence(accs, accs_in_hiconf_space))
        # self.pacbio_hybrid_peptide_set = set(self.pacbio_hybrid_peptide['Base Sequence'].unique())
#***************************************************
# Metamorpheus Protein Groups
#***************************************************
class ProteinGroup:
    def __init__(self) -> None:
        self.gencode_protein_group = read_protein_group_file(gencode_protein_group_file, gencode_gene_map)
        self.gencode_protein_group['is_high_confidence'] = self.gencode_protein_group['genes'].apply(lambda genes: is_high_confidence(genes, genes_in_hiconf_space))

        self.uniprot_protein_group = read_protein_group_file(uniprot_protein_group_file, uniprot_gene_map)
        self.uniprot_protein_group['is_high_confidence'] = self.uniprot_protein_group['genes'].apply(lambda genes: is_high_confidence(genes, genes_in_hiconf_space))

        # pacbio_filtered_protein_group = read_protein_group_file(pacbio_filtered_protein_group_file, pacbio_gene_map)
        self.pacbio_hybrid_protein_group = read_protein_group_file(pacbio_hybrid_protein_group_file, hybrid_gene_map)
        self.pacbio_hybrid_protein_group['is_high_confidence'] = self.pacbio_hybrid_protein_group['accs'].apply(lambda accs: is_high_confidence(accs, accs_in_hiconf_space))

        self.pacbio_rescue_resolve_protein_group = read_protein_group_file(pacbio_rescue_resolve_protein_group_file, hybrid_gene_map)
        self.pacbio_rescue_resolve_protein_group['is_high_confidence'] = self.pacbio_rescue_resolve_protein_group['accs'].apply(lambda accs: is_high_confidence(accs, accs_in_hiconf_space))


