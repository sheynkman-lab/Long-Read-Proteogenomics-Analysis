# get number of gencode transcripts in comprehensive vs basic annotations


# %%
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import config
gtf_file = f'{config.REFERENCE_DIRECTORY}/gencode.v35.annotation.gtf'

# read in gencode gtf

def is_protein_coding(acc_line):
    if 'transcript_type "protein_coding"' in acc_line:
        return True

def is_basic_transcript(acc_line):
    if 'tag "basic"' in acc_line:
        return '1'
    else:
        return '0'

def get_isoname(acc_line):
    return acc_line.split('transcript_name "')[1].split('"')[0]

def get_gene(acc_line):
    return acc_line.split('gene_name "')[1].split('"')[0]


info = [] # list of [gene, isoname, is_basic]
for line in open(gtf_file):
    if line.startswith('#'): continue
    wds = line.split('\t')
    feat = wds[2]
    if feat == 'transcript':
        acc_line = wds[8]
        if is_protein_coding(acc_line):
            gene = get_gene(acc_line)
            isoname = get_isoname(acc_line)
            is_basic = is_basic_transcript(acc_line)
            info.append([gene, isoname, is_basic])

with open('gene_isoname_isbasic.tsv', 'w') as ofile:
    ofile.write('gene\tisoname\tis_basic\n')
    for line in info:
        ofile.write('\t'.join(line) + '\n')



