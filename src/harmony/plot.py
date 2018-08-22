import warnings
import os
import numpy as np
import pandas as pd
import copy
from itertools import chain
from sklearn.preprocessing import StandardScaler

import matplotlib
from matplotlib import font_manager

try:
    os.environ['DISPLAY']
except KeyError:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

with warnings.catch_warnings():
    # catch experimental ipython widget warning
    warnings.simplefilter('ignore')
    import seaborn as sns
    sns.set(context="paper", style='ticks',
            font_scale=1.5, font='Bitstream Vera Sans')

# set plotting defaults
with warnings.catch_warnings():
    # catch warnings that system can't find fonts
    warnings.simplefilter('ignore')
    import seaborn as sns
    fm = font_manager.fontManager
    fm.findfont('Raleway')
    fm.findfont('Lato')

warnings.filterwarnings(action="ignore", message="remove_na is deprecated")

from fa2 import ForceAtlas2
import pandas as pd
import numpy as np
import random

def force_directed_layout(affinity_matrix, cell_names=None, verbose=True, iterations=500):
    """" Function to compute force directed layout from the affinity_matrix
    :param affinity_matrix: Sparse matrix representing affinities between cells
    :param cell_names: pandas Series object with cell names
    :param verbose: Verbosity for force directed layout computation
    :param iterations: Number of iterations used by ForceAtlas 
    :return: Pandas data frame representing the force directed layout
    """

    init_coords = np.random.random((affinity_matrix.shape[0], 2))

    forceatlas2 = ForceAtlas2(
        # Behavior alternatives
        outboundAttractionDistribution=False,  
        linLogMode=False,  
        adjustSizes=False,  
        edgeWeightInfluence=1.0,
        # Performance
        jitterTolerance=1.0,  
        barnesHutOptimize=True,
        barnesHutTheta=1.2,
        multiThreaded=False,  
        # Tuning
        scalingRatio=2.0,
        strongGravityMode=False,
        gravity=1.0,
        # Log
        verbose=verbose)

    positions = forceatlas2.forceatlas2(
        affinity_matrix, pos=init_coords, iterations=iterations)
    positions = np.array(positions)

    # Convert to dataframe
    if cell_names is None:
        cell_names = np.arange(affinity_matrix.shape[0])

    positions = pd.DataFrame(positions,
                             index=cell_names, columns=['x', 'y'])
    return positions


class FigureGrid:
    """
    Generates a grid of axes for plotting
    axes can be iterated over or selected by number. e.g.:
    >>> # iterate over axes and plot some nonsense
    >>> fig = FigureGrid(4, max_cols=2)
    >>> for i, ax in enumerate(fig):
    >>>     plt.plot(np.arange(10) * i)
    >>> # select axis using indexing
    >>> ax3 = fig[3]
    >>> ax3.set_title("I'm axis 3")
    """

    # Figure Grid is favorable for displaying multiple graphs side by side.

    def __init__(self, n: int, max_cols=3, scale=3):
        """
        :param n: number of axes to generate
        :param max_cols: maximum number of axes in a given row
        """

        self.n = n
        self.nrows = int(np.ceil(n / max_cols))
        self.ncols = int(min((max_cols, n)))
        figsize = self.ncols * scale, self.nrows * scale

        # create figure
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=self.ncols)
        self.figure = plt.figure(figsize=figsize)

        # create axes
        self.axes = {}
        for i in range(n):
            row = int(i // self.ncols)
            col = int(i % self.ncols)
            self.axes[i] = plt.subplot(self.gs[row, col])

    def __getitem__(self, item):
        return self.axes[item]

    def __iter__(self):
        for i in range(self.n):
            yield self[i]


def plot_timepoints(layout, timepoints):
    """Plot timepoints on the force directed layoug
    :param layout: Force directed layout
    :param timepoints: Pandas series of timepoints
    """

    # Cluster colors
    n_clusters = len(set(timepoints))
    cluster_colors = pd.Series(sns.color_palette(
        'muted', n_clusters).as_hex(), index=np.sort(timepoints.unique()))

    # Set up figure
    fig = FigureGrid(n_clusters, 5)
    for t, ax in zip(cluster_colors.index, fig):
        ax.scatter(layout['x'], layout['y'], s=3, color='lightgrey')
        cells = timepoints.index[timepoints == t]
        ax.scatter(layout.loc[cells, 'x'], layout.loc[cells, 'y'], 
            s=5, color=cluster_colors[t])
        ax.set_axis_off()
        ax.set_title(t)


def plot_tp_gene_expression(data, layout, genes, timepoints):
    """ Plot gene expression on force directed layout
    :param genes: Iterable of strings to plot on force directed layout
    """

    not_in_dataframe = set(genes).difference(data.columns)
    if not_in_dataframe:
        if len(not_in_dataframe) < len(genes):
            print('The following genes were either not observed in the experiment, '
                  'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
        else:
            print('None of the listed genes were observed in the experiment, or the '
                  'wrong symbols were used.')
            return

    # remove genes missing from experiment
    genes = set(genes).difference(not_in_dataframe)

    # Plot
    ncols = len(genes); nrows = len(timepoints.unique());
    fig = plt.figure(figsize=[4 * ncols, 4 * nrows])
    gs = plt.GridSpec(nrows, ncols)

    for i, t in enumerate(np.sort(timepoints.unique())):
        cells = timepoints.index[timepoints == t]

        for j, gene in enumerate(genes):
            c = data.loc[:, gene]

            ax = plt.subplot(gs[i, j])
            ax.scatter(layout.loc[:, 'x'], layout.loc[:, 'y'], s=3, color='lightgrey')
            ax.scatter(layout.loc[cells, 'x'], layout.loc[cells, 'y'], s=5,
                cmap=matplotlib.cm.Spectral_r, c=c[cells], vmin=c.min(), vmax=c.max())
            ax.set_title(f'{t}: {gene}')
            ax.set_axis_off()


