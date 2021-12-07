#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as tick
from matplotlib.ticker import ScalarFormatter
# matplotlib font parameters
font = {'family' : 'sans-serif',
        'sans-serif':['Arial'],
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)

def reformat_large_tick_values(tick_val, pos):
    """
    Turns large tick values (in the billions, millions and thousands) such as 4500 into 4.5K and also appropriately turns 4000 into 4K (no zero after the decimal).
    """
    if tick_val >= 1000000000:
        val = round(tick_val/1000000000, 1)
        new_tick_format = '{:}B'.format(val)
    elif tick_val >= 1000000:
        val = round(tick_val/1000000, 1)
        new_tick_format = '{:}M'.format(val)
    elif tick_val >= 1000:
        val = round(tick_val/1000, 1)
        new_tick_format = '{:}K'.format(val)
    elif tick_val < 1000:
        new_tick_format = round(tick_val, 1)
    else:
        new_tick_format = tick_val

    # make new_tick_format into a string value
    new_tick_format = str(new_tick_format)
    
    # code below will keep 4.5M as is but change values such as 4.0M to 4M since that zero after the decimal isn't needed
    index_of_decimal = new_tick_format.find(".")
    
    if index_of_decimal != -1:
        value_after_decimal = new_tick_format[index_of_decimal+1]
        if value_after_decimal == "0":
            # remove the 0 after the decimal point since it's not needed
            new_tick_format = new_tick_format[0:index_of_decimal] + new_tick_format[index_of_decimal+2:]
            
    return new_tick_format

def subsample_category_plot(
    subsample, 
    ofile, 
    categories_to_keep=None, 
    ylabel='Number of detected genes', 
    xlabel='Number of subsampled reads', 
    title='Rarefaction, by category',
    normalize = False,
    scale = 'standard'):
    
    if categories_to_keep is not None:
        subsample = subsample[subsample['category'].isin(categories_to_keep)]
    fig, ax = plt.subplots(figsize=(8,5))
    if normalize:
        subsample['mean'] = subsample['mean'] / subsample.groupby('category')['mean'].transform(np.max)

    
    sns.lineplot(data=subsample, y = 'mean', x='size', hue='category', ax=ax)
    if scale == 'log':
        ax.set(yscale='log')

    ax.xaxis.set_major_formatter(tick.FuncFormatter(reformat_large_tick_values))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.title(title)
    plt.savefig(ofile, bbox_inches='tight')


def subsample_plot(gene_ifile, isoform_ifile, ofile, normalize=False, scale='standard'):
    gene_subsample = pd.read_csv(gene_ifile, delimiter=' ', skiprows=1 )
    isoform_subsample = pd.read_csv(isoform_ifile, delimiter=' ', skiprows=1 )
    fig, ax = plt.subplots(figsize=(8,5))

    if normalize:
        gene_subsample['mean'] = gene_subsample['mean'] / gene_subsample['mean'].max()
        isoform_subsample['mean'] = isoform_subsample['mean'] / isoform_subsample['mean'].max()

    sns.lineplot(data=gene_subsample, y = 'mean', x='size', ax=ax, label='gene')
    sns.lineplot(data=isoform_subsample, y = 'mean', x='size', ax=ax, label='isoform')
    # ax.ticklabel_format(style='plain') 
    ax.xaxis.set_major_formatter(tick.FuncFormatter(reformat_large_tick_values))
    if scale == 'log':
        ax.set(yscale='log')
    if normalize:
        ax.set_ylabel(f'Fraction of detected genes/isoforms')
    else:
        ax.set_ylabel(f'Number of detected genes/isoforms')

    # plt.xticks(rotation=45)
    ax.set_xlabel('Number of subsampled reads')
    
    plt.title('Rarefaction, gene or isoform level')
    plt.savefig(ofile, bbox_inches='tight')


# %%
for normalize in [True, False]:
    for scale in ['standard', 'log']:
        ylabel = 'Fraction' if normalize else 'Number'
        isoform_ifile = './jurkat.rarefaction.by_pbid.min_fl_2.step_50k.txt'
        gene_ifile = './jurkat.rarefaction.by_refgene.min_fl_2.step_50k.txt'
        subsample_plot(gene_ifile, isoform_ifile, f'./plots/jurkat.rarefaction.gene_isoform_rarefaction_plot.normalize_{normalize}.scale_{scale}.pdf', normalize, scale)

        isoform_length_subsample = pd.read_csv('./jurkat.rarefaction.by_pbid.min_fl_2.step_50k.length_category.txt', sep=' ', skiprows=1)
        subsample_category_plot(
            isoform_length_subsample,
            f'plots/jurkat.rarefaction.by_pbid.min_fl_2.step_50k.length_category.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=None,
            ylabel=f'{ylabel} of isoforms detected',
            title='Rarefaction, by length category',
            normalize=normalize,
            scale=scale)

        gene_length_subample = pd.read_csv('./jurkat.rarefaction.by_refgene.min_fl_2.step_50k.length_category.txt', sep=' ', skiprows=1)
        subsample_category_plot(
            gene_length_subample,
            f'plots/jurkat.rarefaction.by_refgene.min_fl_2.step_50k.length_category.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=None,
            ylabel=f'{ylabel} of genes detected',
            title='Rarefaction, by length category',
            normalize=normalize,
            scale=scale)

        isoform_length_subsample_exp = pd.read_csv('./jurkat.rarefaction.by_pbid.min_fl_2.step_50k.length_category_expanded.txt', sep=' ', skiprows=1)
        subsample_category_plot(
            isoform_length_subsample_exp,
            f'plots/jurkat.rarefaction.by_pbid.min_fl_2.step_50k.length_category_expanded.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=None,
            ylabel=f'{ylabel} of isoforms detected',
            title='Rarefaction, by length category',
            normalize=normalize,
            scale=scale)

        gene_length_subsample_exp = pd.read_csv('./jurkat.rarefaction.by_refgene.min_fl_2.step_50k.length_category_expanded.txt', sep=' ', skiprows=1)
        subsample_category_plot(
            gene_length_subsample_exp,
            f'plots/jurkat.rarefaction.by_refgene.min_fl_2.step_50k.length_category_expanded.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=None,
            ylabel=f'{ylabel} of genes detected',
            title='Rarefaction, by length category',
            normalize=normalize,
            scale=scale)


        isoform_sqanti_subsample = pd.read_csv('./jurkat.rarefaction.by_pbid.min_fl_2.step_50k.sqanti_category.txt', sep=' ', skiprows=1)
        subsample_category_plot(
            isoform_sqanti_subsample,
            f'plots/jurkat.rarefaction.by_pbid.min_fl_2.step_50k.sqanti_category.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=['full-splice_match', 'incomplete-splice_match', 'novel_not_in_catalog', 'novel_in_catalog'],
            ylabel=f'{ylabel} of isoforms detected',
            title='Rarefaction, by SQANTI category',
            normalize=normalize,
            scale=scale)

        gene_sqanti_subsample = pd.read_csv('./jurkat.rarefaction.by_refgene.min_fl_2.step_50k.sqanti_category.txt', sep=' ', skiprows=1)
        subsample_category_plot(
            gene_sqanti_subsample,
            f'plots/jurkat.rarefaction.by_refgene.min_fl_2.step_50k.sqanti_category.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=['full-splice_match', 'incomplete-splice_match', 'novel_not_in_catalog', 'novel_in_catalog'],
            ylabel=f'{ylabel} of genes detected',
            title='Rarefaction, by SQANTI category',
            normalize=normalize,
            scale=scale)

        peptide_subsample= pd.read_table('./jurkat.filtered.peptides_subsampled.step_25000.iter_20.tsv')
        peptide_subsample['category'] = 'peptides'
        subsample_category_plot(
            peptide_subsample,
            f'plots/jurkat.filtered.peptides_subsampled.step_25000.iter_20.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=None,
            ylabel=f'{ylabel} of peptides detected in mass-spec',
            title='Rarefaction, peptide detection',
            normalize=normalize,
            scale=scale)

        isoform_length_subsample = pd.read_csv('./jurkat.rarefaction.by_pbid.min_fl_2.step50k.length_category_to_10kb.txt', sep=' ', skiprows=1)
        subsample_category_plot(
            isoform_length_subsample,
            f'plots/jurkat.rarefaction.by_pbid.min_fl_2.step_50k.length_category_10kb.normalize_{normalize}.scale_{scale}.pdf',
            categories_to_keep=None,
            ylabel=f'{ylabel} of isoforms detected',
            title='Rarefaction, by length category',
            normalize=normalize,
            scale=scale)
    # %%
