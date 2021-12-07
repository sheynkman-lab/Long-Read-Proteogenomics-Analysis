python make_file_for_subsampling_from_collapsed.py -i isoseq/jurkat.collapsed -o isoseq/jurkat.pre_subsample -m2 isoseq/jurkat_classification.txt
#*********************************************
# LENGTH
#*********************************************
# pbid
python subsample_with_category.py --by pbid --min_fl_count 2 --step 50000 \
subset_data/jurkat.length_categorized.txt > jurkat.rarefaction.by_pbid.min_fl_2.step_50k.length_category.txt

python subsample_with_category.py --by pbid --min_fl_count 2 --step 50000 \
subset_data/jurkat.length_categorized_expanded.txt > jurkat.rarefaction.by_pbid.min_fl_2.step_50k.length_category_expanded.txt
# refgene
python subsample_with_category.py --by refgene --min_fl_count 2 --step 50000 \
subset_data/jurkat.length_categorized.txt > jurkat.rarefaction.by_refgene.min_fl_2.step_50k.length_category.txt

python subsample_with_category.py --by refgene --min_fl_count 2 --step 50000 \
subset_data/jurkat.length_categorized_expanded.txt > jurkat.rarefaction.by_refgene.min_fl_2.step_50k.length_category_expanded.txt

#*********************************************
# SQANTI CATEGORY
#*********************************************
# pbid
python subsample_with_category.py --by pbid --min_fl_count 2 --step 50000 \
subset_data/jurkat.sqanti_category.txt > jurkat.rarefaction.by_pbid.min_fl_2.step_50k.sqanti_category.txt
# refgene
python subsample_with_category.py --by refgene --min_fl_count 2 --step 50000 \
subset_data/jurkat.sqanti_category.txt > jurkat.rarefaction.by_refgene.min_fl_2.step_50k.sqanti_category.txt
#*********************************************
# NO CATEGORY
#*********************************************
# pbid
python subsample.py --by pbid --min_fl_count 2 --step 50000 \
subset_data/jurkat.sqanti_category.txt > jurkat.rarefaction.by_pbid.min_fl_2.step_50k.txt 
# refgene
python subsample.py --by refgene --min_fl_count 2 --step 50000 \
subset_data/jurkat.sqanti_category.txt > jurkat.rarefaction.by_refgene.min_fl_2.step_50k.txt 

