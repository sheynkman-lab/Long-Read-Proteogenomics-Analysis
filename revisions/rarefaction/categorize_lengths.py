#%%
import pandas as pd
import numpy as np

def set_length_category(length, length_categories):
    for category, length_range in length_categories.items():
        if length_range[0] <= length < length_range[1]:
            return category
    return 'NA'


def categorize_lengths(ifile, ofile, length_categories):
    sample = pd.read_table(ifile)
    sample['category'] = sample['length'].apply(lambda length: set_length_category(length, length_categories))
    sample.to_csv(ofile, sep='\t', index=False)

# %%

ifile = './subset_data/jurkat.sqanti_category.txt'

# ofile = './subset_data/jurkat.length_categorized.txt'
# length_categories = {
#     '0-1kb': [0, 1000],
#     '1kb-4kb': [1000, 4000],
#     '4kb+': [4000, np.inf]
# }
# categorize_lengths(ifile, ofile, length_categories)
# length_categories = {
#     '0-1kb': [0, 1000],
#     '1kb-2kb': [1000, 2000],
#     '2kb-3kb': [2000, 3000],
#     '3kb-4kb': [3000, 4000],
#     '4kb-5kb': [4000, 5000],
#     '5kb+': [5000, np.inf]
# }

length_categories = {
    '0-1kb': [0, 1000],
    '1kb-2kb': [1000, 2000],
    '2kb-3kb': [2000, 3000],
    '3kb-4kb': [3000, 4000],
    '4kb-5kb': [4000, 5000],
    '5kb-6kb': [5000, 6000],
    '6kb-7kb': [6000, 7000],
    '7kb-8kb': [7000, 8000],
    '8kb-9kb': [8000, 9000],
    '9kb-10kb': [9000, 10000],
    '10kb+': [10000, np.inf]
}
ofile = './subset_data/jurkat.length_categorized_10kb.txt'

categorize_lengths(ifile, ofile, length_categories)
# %%
