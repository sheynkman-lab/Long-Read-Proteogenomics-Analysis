
GMST_PROG=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/SQANTI3/utilities/gmst/gmst.pl
ODIR=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics-Analysis/revisions/gmst/results
INPUT_FASTA=/Users/bj8th/Documents/Sheynkman-Lab/Data/JURKAT_06-06-2021/sqanti3-filtered/jurkat_corrected.5degfilter.fasta

perl $GMST_CMD --faa --strand direct --fnn --output $ODIR $INPUT_FASTA


 