import matplotlib.pyplot as plt
import numpy as np

def mask_first(x):
    result = np.ones_like(x)
    result[0] = 0
    return result

    
def get_major_isoform(sqanti_info):
    major_isoform = sqanti_info \
        .sort_values(by = 'cpm', ascending = False) \
        .groupby(['gene']).first()
    return major_isoform

def get_minor_isoforms(sqanti_info):
    mask = sqanti_info \
        .sort_values(by = 'cpm', ascending = False) \
        .groupby(['gene'])['gene'] \
            .transform(mask_first).astype(bool)
    minor_isoforms = sqanti_info.loc[mask]
    return minor_isoforms

def plot_appris_donut(major_isoform, figure):
    # donut chart, how many major isoforms are appris versus non-appris
    # read in appris isonames (note that some genes have two isonames as best appris, both are listed here)
    
    with open(f'stats/{figure}_number_major_isoforms_that_match_appris.txt', 'w') as ofile:
        major_isoform_is_appris = major_isoform[major_isoform['is_appris']]
        major_isoform_is_not_appris = major_isoform[~major_isoform['is_appris']]
        num_is_appris = major_isoform_is_appris.shape[0]
        num_is_not_appris = major_isoform_is_not_appris.shape[0]
        ofile.write(f'Number of major isoforms is appris:    {num_is_appris} ({num_is_appris/(num_is_appris+num_is_not_appris)*100:.1f}%)\n')
        ofile.write(f'Number of major isoforms is not appris:{num_is_not_appris} ({num_is_not_appris/(num_is_appris+num_is_not_appris)*100:.1f}%)')
    # appris pie / donut chart
    appris_labels = ['major isoform\nconsidered reference', 'major isoform\nnot considered reference']
    appris_sizes = [
        num_is_appris,
        num_is_not_appris,
    ]
    fig, ax = plt.subplots(figsize=(6,6))
    colors = ['#66b3ff','#ff9999']
    def label_donut(pct, allvals):
        absolute = round(pct/100.*np.sum(allvals))
        return "{:.1f}%\n({:,})".format(pct, absolute)
    wedges, texts, autotexts = ax.pie(appris_sizes, autopct=lambda pct: label_donut(pct, appris_sizes), colors=colors,startangle=90)
    centre_circle = plt.Circle((0,0),0.40,fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)# Equal aspect ratio ensures that pie is drawn as a circle
    ax.axis('equal')  
    ax.legend(wedges, ['Yes','No'], 
                title='Major isoform is considered\nthe reference (APPRIS)',
                loc="center left",
                bbox_to_anchor=(0.25, 0, 0.5, -0.4),
                fontsize=20)
    plt.tight_layout()
    plt.savefig(f'plot/{figure}_donut_number_of_major_appris.pdf', bbox_inches='tight')
    plt.show()
    plt.clf()