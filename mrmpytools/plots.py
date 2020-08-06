"""
plots.py
Some functions for plotting
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import string ## used for panel labeling (A, B, C, ...)

def sci_fmt(val):
    """
    Format a value in scientific notation.
    Intendend for manually formatting the tick labels
    on axes of plots.
    
    Examples: 
        0.2 -> $2 \cdot 10^{-1}$
        300 -> $3 \cdot 10^{2}$
    """
    exp = np.floor(np.log10(val))
    mul = np.floor(val / 10**exp)
    sci_str = f"${mul:0.0f} \\cdot 10^{{{exp:0.0f}}}$"
    return sci_str


def ordinate_transform(ax, z, which='x'):
    """in order to get coordinates in a subplot,
    one can use transforms. This function does transforms
    for only one ordinate.
    TODO: this is obsolete, use transform factory
    """
    inv = ax.transData.inverted()
    if which == 'x':
        u, v = ax.transAxes.transform((z, 0))
        x, y = inv.transform((u, 0))
        return x
    elif which == 'y':
        u, v = ax.transAxes.transform((0, z))
        x, y = inv.transform((0, v))
        return y
    else:
        raise Exception("which must be either 'x' or 'y'")

def color_violins(vp, color, alpha=1):
    """
    auxilliary function to color the bodies of violinplots
    """
    for body in vp["bodies"]:
        body.set_alpha(alpha)
        body.set_color(color)
        
def violin_matrix_plot(subjects, parnames, mu_parnames, pretty_parnames, pardims, chain, **kwargs):
    """
    plot individual parameter estimates and posterior predictive distributions
    in a single comprehensive figure.
    FIXME: replace hyper with posterior predictive distribution
    """
    NumSubjects = len(subjects)
    gs = GridSpec(len(parnames), NumSubjects+1)
    gs.update(hspace=0.25)

    fig = plt.figure(figsize=(7,2.5*len(parnames)))
    axs, bxs = [], []

    for i, parname in enumerate(parnames):
        ax = fig.add_subplot(gs[i,:NumSubjects])
        axs.append(ax)
        data = chain[parname]
        if pardims[i]==1:
            violins = ax.violinplot(data, positions=range(NumSubjects), **kwargs)
            color_violins(violins, 'tab:blue', alpha=0.5)
        else:
            twined_data = np.concatenate([data[:,:,j] for j in range(pardims[i])], axis=1)
            pos = np.array(range(NumSubjects))
            pos = np.concatenate([pos + dp for dp in np.linspace(-0.25, 0.25, pardims[i])])
            violins = ax.violinplot(twined_data, positions=pos, **kwargs)
            color_violins(violins, 'tab:blue', alpha=0.5)
            ## TODO: give violins different colors
        ax.set_ylabel(pretty_parnames[i])
        ax.set_xticks(range(NumSubjects))
        ax.set_xticklabels(subjects)
        if mu_parnames[i] is not None:
            bx = fig.add_subplot(gs[i, NumSubjects], sharey=ax)
            data = chain[mu_parnames[i]]
            if pardims[i]==1:
                violins = bx.violinplot(data, **kwargs)
                color_violins(violins, 'tab:blue', alpha=0.5)
            else:
                pos = 1 + np.linspace(-0.25, 0.25, pardims[i])
                violins = bx.violinplot(data, positions=pos, **kwargs)
                color_violins(violins, 'tab:blue', alpha=0.5)
            bx.set_xticks([1])
            bx.set_xticklabels(["hyper"])
            bx.axes.get_yaxis().set_visible(False)
            bxs.append(bx)
        else:
            bxs.append(None)
        ## label for panel
        if len(parnames) > 1:
            ax.text(-0.15, 0.9, string.ascii_uppercase[i], fontsize='x-large', 
                    ha='right', va='bottom',transform=ax.transAxes)

    axs[0].get_shared_x_axes().join(*axs)
    axs[0].autoscale(axis='x')
    return fig, axs, bxs


def simple_boxplot(ax, pos, data, color='k', color_med='red', p=50, horizontal=False, s=None, **kwargs):
    uppers = [np.nanpercentile(x, 50+p/2) for x in data]
    downers = [np.nanpercentile(x, 50-p/2) for x in data]
    medians = [np.nanmedian(x) for x in data]
    if horizontal:
        ax.scatter(uppers, pos, marker='|', color=color, s=s, **kwargs)
        ax.scatter(downers, pos, marker='|', color=color, s=s, **kwargs)
        ax.scatter(medians, pos, marker='|', color=color_med, s=s, **kwargs)
        for p, d, u in zip(pos, downers, uppers):
            ax.plot([d, u], [p, p], color=color, **kwargs)        
    else: ## vertical
        ax.scatter(pos, uppers, marker='_', color=color, s=s, **kwargs)
        ax.scatter(pos, downers, marker='_', color=color, s=s, **kwargs)
        ax.scatter(pos, medians, marker='_', color=color_med, s=s, **kwargs)
        for p, d, u in zip(pos, downers, uppers):
            ax.plot([p, p], [d, u], color=color, **kwargs)
