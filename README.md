# Long-Read-Proteogenomics-Analysis

### Download Reference Data

The reference data can be downloaded [here](https://zenodo.org/record/5076056#.YOSPAhNKhTY)

### Download Results Data
 
To use the data used in the "Enhanced protein isoform characterization through long-read proteogenomics" manuscript download the data from this :

s3://sheynkman-lab-lifebit/deploit/users/5fa6c3c061aba10112169d7d/projects/602eb362e9e4e10112f0dd63/jobs/60bcb29b303ee601a69d8c74/

### Update Config

In the config.py file, update the PIPELINE_RESULTS_DIRECTORY and REFERENCE_DIRECTORY to the paths of the two directories where the data above is downloaded

```

PIPELINE_RESULTS_DIRECTORY='path/to/results/directory'
REFERENCE_DIRECTORY='path/to/reference/directory'
```

### Run Jupyter Notebooks 

Within the 0_pre_analysis directory, run get_appris_reference_isoforms.py and get_number_gencode_transcripts.py

```
python get_appris_reference_isoforms.py  
python get_number_gencode_transcripts.py
```

Open each .ipynb file and run each notebook fully. Figures are stored in the /plot directory and statistics are stored in the /stats directory
