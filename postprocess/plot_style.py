"""Shared style for the editable notebook figures."""

import matplotlib as mpl


def setup():
    mpl.rcParams.update({
        'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsmath}',
        'font.family': 'serif', 'font.size': 20, 'axes.labelsize': 20,
        'axes.titlesize': 20, 'xtick.labelsize': 20, 'ytick.labelsize': 20,
        'legend.fontsize': 20, 'legend.frameon': False, 'axes.grid': False,
        'xtick.direction': 'out', 'ytick.direction': 'out',
    })
