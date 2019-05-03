import pandas as pd
import numpy as np
import os
import re
from itertools import chain
from collections import OrderedDict
import warnings
from sklearn.decomposition import PCA


def load_from_csvs(csv_files, sample_names=None, min_cell_count=10):
    """ Read in scRNA-seq data from csv files
    :param csv_files: List of csv files
    :param sample_names: Prefix to attach to the cell names
    :param min_cell_count: Minimum number of cells a gene should be detected in
    :return: Pandas data frame representing the count matrix
    """
    if sample_names is None:
        sample_names = [re.sub('.csv', '', os.path.split(i)[-1]) for i in csv_files]


    # Load counts
    print('Loading count matrices...')
    counts_dict = OrderedDict()

    for csv_file, sample in zip(csv_files, sample_names):
        print(sample)
        counts = pd.read_csv(csv_file, index_col=0)

        # Update sample names
        counts.index = sample + '_' + counts.index.astype(str)

        # Remove zero count genes
        counts = counts.loc[:, counts.sum() > 0]
        counts = counts.astype(np.int16)

        # Update dictionary
        counts_dict[sample] = counts

    # Concatenate cells and genes
    print('Concatenating data..')
    all_cells = list(chain(*[list(counts_dict[sample].index)
                             for sample in sample_names]))
    all_genes = list(
        chain(*[list(counts_dict[sample].columns) for sample in sample_names]))
    all_genes = list(set(all_genes))

    # Counts matrix
    counts = pd.DataFrame(0, index=all_cells,
                          columns=all_genes, dtype=np.int16)
    for sample in sample_names:
        sample_counts = counts_dict[sample]
        counts.loc[sample_counts.index, sample_counts.columns] = sample_counts

    # Filter out low detection genes
    gs = counts.sum()
    counts = counts.loc[:, counts.columns[gs > min_cell_count]]

    return counts


def hvg_genes(norm_df, no_genes=1000):
    """ Select genes based on normalized dispersion
    """

    # Mean and variance
    mean = norm_df.mean()
    var = norm_df.var()

    # Construct data frame
    df = pd.DataFrame()
    df['mean'] = mean
    df['dispersion'] = var/mean
    df[df.isnull()] = 0

    # Normalized disperion
    from statsmodels import robust
    df['mean_bin'] = pd.cut(df['mean'], np.r_[-np.inf,
                                              np.percentile(df['mean'],
                                                            np.arange(10, 105, 5)), np.inf])
    disp_grouped = df.groupby('mean_bin')['dispersion']
    disp_median_bin = disp_grouped.median()
    # the next line raises the warning: "Mean of empty slice"
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        disp_mad_bin = disp_grouped.apply(robust.mad)
    df['dispersion_norm'] = np.abs((df['dispersion'].values - disp_median_bin[df['mean_bin'].values].values)) \
        / disp_mad_bin[df['mean_bin'].values].values

    # Subset of genes
    use_genes = df['dispersion_norm'].sort_values().index[::-1][:no_genes]

    return use_genes


def run_pca(data, n_components=300, var_explained=0.85):
    """Run PCA

    :param data: Dataframe of cells X genes. Typicaly multiscale space diffusion components
    :param n_components: Number of principal components
    :param var_explained: Include components that explain amount variance. Note
    number of components = min(n_components, components explaining var_explained)
    :return: PCA projections of the data and the explained variance
    """
    init_components = min([n_components, data.shape[0]])
    pca = PCA(n_components=init_components, svd_solver='randomized')
    pca.fit(data)
    if pca.explained_variance_ratio_.sum() >= 0.85:
        n_components = np.where(np.cumsum(pca.explained_variance_ratio_) >= var_explained)[0][0]

    print(f'Running PCA with {n_components} components')
    pca = PCA(n_components=n_components, svd_solver='randomized')
    pca_projections = pca.fit_transform(data)
    pca_projections = pd.DataFrame(pca_projections, index=data.index)
    return pca_projections, pca.explained_variance_ratio_


def normalize_counts(data, multiplier=10000):
    """Correct the counts for molecule count variability

    :param data: Counts matrix: Cells x Genes
    :return: Normalized matrix
    """
    ms = data.sum(axis=1)
    norm_df = data.div(ms, axis=0) * multiplier
    return norm_df


def log_transform(data, pseudo_count=0.1):
    """Log transform the matrix

    :param data: Counts matrix: Cells x Genes
    :return: Log transformed matrix
    """
    return np.log2(data + pseudo_count)
