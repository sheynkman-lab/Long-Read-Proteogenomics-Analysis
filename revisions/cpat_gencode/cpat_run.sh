hexamer=Human_Hexamer.tsv
logit_model=Human_logitModel.RData
sample_fasta=gencode.v35.pc_transcripts.fa
name=GENCODE

cpat.py \
-x $hexamer \
-d $logit_model \
-g $sample_fasta \
--min-orf=50 \
--top-orf=50 \
-o $name \
1> ${name}_cpat.output \
2> ${name}_cpat.error