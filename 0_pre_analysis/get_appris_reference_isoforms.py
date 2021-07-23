# get list of genes and the isonames that represent the appris reference

# %%

from collections import defaultdict
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import config
gtf_file = f'{config.REFERENCE_DIRECTORY}/gencode.v35.annotation.gtf'


# gene -> [(appris, isoname), (appris, isoname)]
appris = defaultdict(list)
for line in open(gtf_file): 
    if line.startswith('#'): continue
    wds = line.split('\t')
    feat = wds[2]
    if feat == 'transcript':
        acc = wds[8]
        gene = acc.split('gene_name "')[1].split('"')[0]
        isoname = acc.split('transcript_name "')[1].split('"')[0]
        appris_tag = ''
        if 'appris_principal' in acc:
            appris_tag = 'appris_principal_' + acc.split('appris_principal')[1].split('"')[0]
            appris[gene].append([appris_tag, isoname])

def get_the_top_appris_isoname(appris_info):
    best_appris_tag = sorted(appris_info)[0][0]
    isonames_w_best_appris_tag = []
    # sometimes two isonames with same appris status
    # so grab them both
    for appris_tag, isoname in appris_info:
        if appris_tag == best_appris_tag:
            isonames_w_best_appris_tag.append(isoname)
    return best_appris_tag, isonames_w_best_appris_tag

with open('./appris_transcripts.tsv', 'w') as ofile:
    ofile.write('gene\tappris_tag\tisoname\n')
    for gene, appris_info in appris.items():
        best_appris_tag, appris_isonames = get_the_top_appris_isoname(appris_info)
        for isoname in appris_isonames:
            ofile.write('{}\t{}\t{}\n'.format(gene, best_appris_tag, isoname))

# %%
