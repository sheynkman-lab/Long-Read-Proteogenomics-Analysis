input_fasta=/Users/bj8th/Documents/Sheynkman-Lab/Data/JURKAT_06-06-2021/sqanti3-filtered/jurkat_corrected.5degfilter.fasta
output_dir=/Users/bj8th/Documents/Sheynkman-Lab/GitHub/Long-Read-Proteogenomics-Analysis/revisions/transdecoder/RESULTS/transdecoder_results
min_orf_size=50
TransDecoder.LongOrfs -t $input_fasta -m $min_orf_size --output_dir $output_dir
TransDecoder.Predict -t $input_fasta --single_best_only --output_dir $output_dir
