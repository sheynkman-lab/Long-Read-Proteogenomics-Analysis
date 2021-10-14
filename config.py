sqanti_colors = {
    "FSM": "#6BAED6",
    "ISM": "#FC8D59",
    "NIC": "#78C679",
    "NNC": "#EE6A50",
    "Genic Genomic":"#969696",
    "Antisense": "#66C2A4",
    "Fusion":"goldenrod1",
    "Intergenic":"darksalmon",
    "Genic Intron":"#41B6C4",
}
sqanti_protein_colors = {
    'pFSM' : '#1974b5',
    'pNIC' : '#4a943e',
    'pNNC' : '#d94c2e' 
}

database_colors = {
    'GENCODE' : 'darkblue',
    'UniProt' : 'darkgreen',
    'PacBio' : 'darkred', 
    'PacBio Rescue Resolve' : 'red'
}

database_colors_light = {
    'GENCODE' : 'lightsteelblue',
    'UniProt' : 'darkseagreen',
    'PacBio' : 'darksalmon', 
    'PacBio Rescue Resolve' : 'pink'
}
transcript_short_classification = {
    'novel_not_in_catalog': 'NNC', 
    'full-splice_match': 'FSM', 
    'novel_in_catalog' : 'NIC',
    'incomplete-splice_match' : 'ISM', 
    'fusion' : 'fusion',
}

# matplotlib font parameters
font = {'family' : 'sans-serif',
        'sans-serif':['Arial'],
        'weight' : 'normal',
        'size'   : 16}

PIPELINE_RESULTS_DIRECTORY='/Users/bj8th/Documents/Sheynkman-Lab/Data/JURKAT_06-06-2021'
# PIPELINE_RESULTS_DIRECTORY='/Users/gloriasheynkman/Documents/research_drive/projects/long_read_proteogenomics/LRPG-Manuscript/data/results/jurkat'
REFERENCE_DIRECTORY='/Users/bj8th/Documents/Sheynkman-Lab/Data/Reference'
EXPERIMENT_NAME='jurkat'
# REFERENCE_DIRECTORY='/Users/gloriasheynkman/Documents/research_drive/projects/long_read_proteogenomics/LRPG-Manuscript/data/input'