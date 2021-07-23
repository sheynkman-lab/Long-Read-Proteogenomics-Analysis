import pandas as pd
from matplotlib_venn import venn2, venn2_circles
import matplotlib.pyplot as plt
import sys, os
import matplotlib

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import config

matplotlib.rc('font', **config.font)

# customized function to make venn diagrams
def make_venn_diagram(left_set, right_set, left_color, right_color, left_label, right_label, cat,space, filename):
    fig, ax = plt.subplots(figsize=(5,5),tight_layout=True)
    vd = venn2([left_set, right_set], set_labels=None, ax=ax, set_colors=(left_color, right_color), alpha=1)
    c = venn2_circles(subsets=[left_set, right_set], linestyle='solid', linewidth=0.5)
    h,l = [],[]
    total = len(left_set.union(right_set))
    for i in ['10','11','01']:
        num = int(vd.get_label_by_id(i).get_text())
        lbl = vd.get_label_by_id(i)
        lbl.set_text(f'{num}\n{num/total:.0%}')
        x, y = lbl.get_position()
        if i == '10':
            lbl.set_position((x-0.2, y))
        if i == '01':
            lbl.set_position((x+0.2, y))
        h.append(vd.get_patch_by_id(i))
    vd.get_patch_by_id('11').set_color('lightgray')
    l = [left_label, "Both", right_label]
    legend_xpos = 0.92
    if left_label == 'UniProt':
        legend_xpos = 0.86
    ax.legend(handles=h,labels=l,title='Database',loc='best', bbox_to_anchor=(legend_xpos, -0.))
    if space != '':

        plt.title(f"{cat} identified via MS search\n{left_label} vs {right_label}\n{space}")
    else:
        plt.title(f"{cat} identified via MS search\n{left_label} vs {right_label}")

    fig.savefig(f"plot/{filename}.pdf", bbox_inches='tight')
    left_size = len(left_set - right_set)
    right_size = len(right_set - left_set)
    both_size = len(left_set.intersection(right_set))

    with open(f'stats/{filename}_sizes.tsv', 'w') as ofile:
        ofile.write(f'Database\tsize\n')
        ofile.write(f'{left_label}\t{left_size}\n')
        ofile.write(f'Overlap\t{both_size}\n')
        ofile.write(f'{right_label}\t{right_size}\n')

    peptide_df = pd.DataFrame(index = left_set.union(right_set), columns = [left_label,right_label ], data=0, )
    peptide_df.index.name = cat
    peptide_df.loc[left_set, left_label] = 1
    peptide_df.loc[right_set, right_label] = 1
    peptide_df.to_csv(f'stats/{filename}_{cat}.tsv', sep='\t')


# customized function to make venn diagrams, when three sectors known
def make_venn_diagram_three_input(left, right, overlap, left_color, right_color, left_label, right_label, cat, space, figure):
    fig, ax = plt.subplots(figsize=(5,5),tight_layout=True)
    vd = venn2(subsets=[left, right, overlap], set_labels=None, ax=ax, set_colors=(left_color, right_color), alpha=1)
    c = venn2_circles(subsets=[left, right, overlap], linestyle='solid', linewidth=0.5)
    h,l = [],[]
    total = left + right + overlap
    for i in ['10','11','01']:
        num = int(vd.get_label_by_id(i).get_text())
        lbl = vd.get_label_by_id(i)
        lbl.set_text(f'{num}\n{num*1./total:.0%}')
        x, y = lbl.get_position()
        if i == '10':
            cx,cy = vd.get_circle_center(0)
            cx = cx - vd.get_circle_radius(0)
            lbl.set_position((cx-0.2, y))
        if i == '01':
            cx,cy = vd.get_circle_center(1)
            cx = cx + vd.get_circle_radius(1)
            lbl.set_position((cx+0.2, y))
        h.append(vd.get_patch_by_id(i))
    vd.get_patch_by_id('11').set_color('lightgray')
    l = [left_label, "Both", right_label]
    legend_xpos = 0.92
    if cat == 'legend_only':
        ax.legend(handles=h,labels=l,title='Protein database source:',loc='best', ncol=3, bbox_to_anchor=(legend_xpos, -0.))
    else:
        pass
        # no legend
        # ax.legend(handles=h,labels=l,title='Protein database source:',loc='best', bbox_to_anchor=(legend_xpos, -0.))
    if space != "":
        plt.title(f"{cat} identified via MS search\n{left_label} vs {right_label}\n{space}")
    else:
        plt.title(f"{cat} identified via MS search\n{left_label} vs {right_label}")

    fig.savefig(f"plot/{figure}.pdf", bbox_inches='tight')

    with open(f'stats/{figure}_sizes.tsv', 'w') as ofile:
        ofile.write(f'Database\tsize\n')
        ofile.write(f'{left_label}\t{left}\n')
        ofile.write(f'Overlap\t{overlap}\n')
        ofile.write(f'{right_label}\t{right}\n')

