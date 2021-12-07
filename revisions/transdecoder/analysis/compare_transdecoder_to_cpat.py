#%%
import matplotlib_venn
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict

def read_transdecoder_nucleotides_file(transdecoder_nucleotides_file):
    transdecorder_sequences = {}
    for record in SeqIO.parse(transdecoder_nucleotides_file, 'fasta'):
        acc = record.id.split('.')[:-1]
        acc = '.'.join(acc)
        transdecorder_sequences[acc] = record.seq
    return transdecorder_sequences

def read_cpat_nucleotides_file(cpat_best_file, cpat_nucleotides_file):
    cpat_best = pd.read_table(cpat_best_file)
    best_orfs = set(cpat_best['ID'])
    cpat_sequences = {}
    for record in SeqIO.parse(cpat_nucleotides_file,'fasta'):
        if record.id in best_orfs:
            accession = record.id.split('_')[0]
            cpat_sequences[accession] = record.seq
    return cpat_sequences

def compare_sequences(cpat, transdecoder):
    comparisons = []
    accessions = set(cpat.keys()).union(set(transdecoder.keys()))
    for accession in accessions:
        accession_comparison = ''
        if accession not in cpat.keys():
            accession_comparison = 'not_found_in_cpat'
        elif accession not in transdecoder.keys():
            accession_comparison = 'not_found_in_transdecoder'
        elif cpat[accession] == transdecoder[accession]:
            accession_comparison = 'same_orf_called'
        else:
            accession_comparison = 'different_orfs_called'
        comparisons.append([accession, accession_comparison])
    comparisons = pd.DataFrame(comparisons, columns=['PacBio Accession','Comparison'])
    return comparisons
#%%
transdecoder_nucleotides_file = '/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics-Analysis/revisions/transdecoder/analysis/jurkat_corrected.5degfilter.fasta.transdecoder.cds'
cpat_best_file = '/Users/bj8th/Documents/Sheynkman-Lab/Data/JURKAT_06-06-2021/cpat/jurkat.ORF_prob.best.tsv'
cpat_nucleotides_file = '/Users/bj8th/Documents/Sheynkman-Lab/Data/JURKAT_06-06-2021/cpat/jurkat.ORF_seqs.fa'

transdecorder_sequences = read_transdecoder_nucleotides_file(transdecoder_nucleotides_file)
cpat_sequences = read_cpat_nucleotides_file(cpat_best_file, cpat_nucleotides_file)

# %%
comparisions = compare_sequences(cpat_sequences, transdecorder_sequences)
# %%
print(cpat_sequences['PB.7237.20'])
print()
print(transdecorder_sequences['PB.7237.20'])
# %%

def get_all_cpat_orfs(cpat_nucleotides_file):
    cpat_sequences = defaultdict(list)
    for record in SeqIO.parse(cpat_nucleotides_file,'fasta'):
        accession = record.id.split('_')[0]
        cpat_sequences[accession].append(str(record.seq))
    return cpat_sequences

def transdecoder_orf_found_in_cpat(accession, transdecoder, cpat_full):
    transdecoder_seq = str(transdecoder[accession])

    return transdecoder_seq in cpat_full[accession]
# cpat_full = get_all_cpat_orfs(cpat_nucleotides_file)
different_orfs = comparisions[comparisions['Comparison'] =='different_orfs_called']
different_orfs['transdecoder_sequence_found_in_cpat'] = different_orfs['PacBio Accession'].apply(lambda acc: transdecoder_orf_found_in_cpat(acc, transdecorder_sequences, cpat_full))
# %%
comparisions
# %%
