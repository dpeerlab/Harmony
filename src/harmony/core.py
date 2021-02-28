import pandas as pd
import numpy as np
import scanpy as sc

from scipy.sparse import find, csr_matrix

from harmony import utils


def augmented_affinity_matrix(
    data_df,
    timepoints,
    timepoint_connections,
    n_neighbors=30,
    n_jobs=-2,
    pc_components=1000,
    device="cpu",
):
    """Function for max min sampling of waypoints

    :param data_df: Normalized data frame. Data frame should be sorted according to the timepoints
    :param timepoints: Panadas series indicating timepoints for each cell in data_df
    :param timepoint_connections: Links between timepoints
    :param n_neighbors: Number of nearest neighbors for graph construction
    :param n_jobs: Nearest Neighbors will be computed in parallel using n_jobs.
    :param pc_components: Minimum number of principal components to use. Specify `None` to use pre-computed components
    :return: Affinity matrix  augmented to mutually nearest neighbors
    """

    # Timepoints nad data_df should in same order
    timepoints = timepoints[data_df.index]
    cell_order = data_df.index

    # Time point cells and indices
    tp_cells = pd.Series()
    tp_offset = pd.Series()
    offset = 0
    for i in timepoints.unique():
        tp_offset[i] = offset
        tp_cells[i] = list(timepoints.index[timepoints == i])
        offset += len(tp_cells[i])

    # Run PCA to denoise the dropouts
    if pc_components is None:
        pca_projections = data_df
    else:
        pca_projections, _ = utils.run_pca(data_df,device, n_components=pc_components)

    # Nearest neighbor graph construction and affinity matrix
    print('Nearest neighbor computation...')

    # --------------------------------------------------------------------------
    # nbrs = NearestNeighbors(n_neighbors=n_neighbors,
    #                         metric='euclidean', n_jobs=-2)
    # nbrs.fit(pca_projections.values)
    # dists, _ = nbrs.kneighbors(pca_projections.values)
    # adj = nbrs.kneighbors_graph(pca_projections.values, mode='distance')
    # # Scaling factors for affinity matrix construction
    # ka = np.int(n_neighbors / 3)
    # scaling_factors = pd.Series(dists[:, ka], index=cell_order)
    # # Affinity matrix
    # nn_aff = _convert_to_affinity(adj, scaling_factors, True)
    # --------------------------------------------------------------------------
    
    if device == "cpu":
        temp = sc.AnnData(data_df.values)
        sc.pp.neighbors(temp, n_pcs=0, n_neighbors=n_neighbors)
        # maintaining backwards compatibility to Scanpy `sc.pp.neighbors`
        try:
            kNN = temp.uns['neighbors']['distances']
        except KeyError:
            kNN = temp.obsp['distances']
    elif device == "gpu":
        from cuml.neighbors import NearestNeighbors
        nn = NearestNeighbors(n_neighbors=n_neighbors,metric="euclidean")
        X_contiguous = np.ascontiguousarray(data_df.values)
        nn.fit(X_contiguous)

        kNN = nn.kneighbors_graph(X_contiguous,mode="distance")
        kNN.setdiag(0)
        kNN.eliminate_zeros()

    # Adaptive k
    adaptive_k = int(np.floor(n_neighbors / 3))
    scaling_factors = np.zeros(data_df.shape[0])

    for i in np.arange(len(scaling_factors)):
        scaling_factors[i] = np.sort(kNN.data[kNN.indptr[i]:kNN.indptr[i + 1]])[adaptive_k - 1]

    scaling_factors = pd.Series(scaling_factors, index=cell_order)

    # Affinity matrix
    nn_aff = _convert_to_affinity(kNN, scaling_factors, device, True)

    # Mututally nearest neighbor affinity matrix
    # Initilze mnn affinity matrix
    N = len(cell_order)
    full_mnn_aff = csr_matrix(([0], ([0], [0])), [N, N])
    for i in timepoint_connections.index:
        t1, t2 = timepoint_connections.loc[i, :].values
        print(f'Constucting affinities between {t1} and {t2}...')

        # MNN matrix  and distance to ka the distance
        t1_cells = tp_cells[t1]
        t2_cells = tp_cells[t2]
        mnn = _construct_mnn(t1_cells, t2_cells, pca_projections,
                             n_neighbors, device, n_jobs)

        # MNN Scaling factors
        # Distance to the adaptive neighbor
        ka_dists = pd.Series(0.0, index=t1_cells + t2_cells)
        # T1 scaling factors
        ka_dists[t1_cells] = _mnn_ka_distances(mnn, n_neighbors)
        # T2 scaling factors
        ka_dists[t2_cells] = _mnn_ka_distances(mnn.T, n_neighbors)

        # Scaling factors
        mnn_scaling_factors = pd.Series(0.0, index=cell_order)
        mnn_scaling_factors[t1_cells] = _mnn_scaling_factors(
            ka_dists[t1_cells], scaling_factors, device)
        mnn_scaling_factors[t2_cells] = _mnn_scaling_factors(
            ka_dists[t2_cells], scaling_factors, device)

        # MNN affinity matrix
        full_mnn_aff = full_mnn_aff + \
            _mnn_affinity(mnn, mnn_scaling_factors,
                          tp_offset[t1], tp_offset[t2], device)

    # Symmetrize the affinity matrix and return
    aff = nn_aff + nn_aff.T + full_mnn_aff + full_mnn_aff.T
    return aff, nn_aff + nn_aff.T


def _convert_to_affinity(adj, scaling_factors, device, with_self_loops=False):
    """ Convert adjacency matrix to affinity matrix
    """
    N = adj.shape[0]
    rows, cols, dists = find(adj)
    if device == "gpu":
        import cupy as cp
        from cupyx.scipy.sparse import csr_matrix as csr_matrix_gpu
        dists = cp.array(dists) ** 2/(cp.array(scaling_factors.values[rows]) ** 2)
        rows, cols = cp.array(rows), cp.array(cols)
        # Self loops
        if with_self_loops:
            dists = cp.append(dists, cp.zeros(N))
            rows = cp.append(rows, range(N))
            cols = cp.append(cols, range(N))
        aff = csr_matrix_gpu((cp.exp(-dists), (rows, cols)), shape=(N, N)).get()
    elif device == "cpu":
        dists = dists ** 2/(scaling_factors.values[rows] ** 2)

        # Self loops
        if with_self_loops:
            dists = np.append(dists, np.zeros(N))
            rows = np.append(rows, range(N))
            cols = np.append(cols, range(N))
        aff = csr_matrix((np.exp(-dists), (rows, cols)), shape=[N, N])
    return aff


def _construct_mnn(t1_cells, t2_cells, data_df, n_neighbors,device,n_jobs=-2):
    # FUnction to construct mutually nearest neighbors bewteen two points
    
    if device == "gpu":
        from cuml import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=n_neighbors,
                                metric='euclidean')
    elif device == "cpu":
        from sklearn.neighbors import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=n_neighbors,
                                metric='euclidean', n_jobs=n_jobs)
    
    print(f't+1 neighbors of t...')
    nbrs.fit(data_df.loc[t1_cells, :].values)
    t1_nbrs = nbrs.kneighbors_graph(
        data_df.loc[t2_cells, :].values, mode='distance')

    print(f't neighbors of t+1...')
    nbrs.fit(data_df.loc[t2_cells, :].values)
    t2_nbrs = nbrs.kneighbors_graph(
        data_df.loc[t1_cells, :].values, mode='distance')

    # Mututally nearest neighbors
    mnn = t2_nbrs.multiply(t1_nbrs.T)
    mnn = mnn.sqrt()
    return mnn

# Rewritten for speed improvements
def _mnn_ka_distances(mnn, n_neighbors):
    # Function to find distance ka^th neighbor in the mutual nearest neighbor matrix
    ka = np.int(n_neighbors / 3)
    ka_dists = np.repeat(None, mnn.shape[0])
    x, y, z = find(mnn)
    rows=pd.Series(x).value_counts()
    for r in rows.index[rows >= ka]:
        ka_dists[r] = np.sort(z[x==r])[ka - 1]
    return ka_dists

def _mnn_scaling_factors(mnn_ka_dists, scaling_factors,device):
    if device == "gpu":
        from cuml import LinearRegression
    else:
        from sklearn.linear_model import LinearRegression
    cells = mnn_ka_dists.index[~mnn_ka_dists.isnull()]

    # Linear model fit
    x = scaling_factors[cells]
    y = mnn_ka_dists[cells]
    lm = LinearRegression()
    lm.fit(x.values.reshape(-1, 1), y.values.reshape(-1, 1))

    # Predict
    x = scaling_factors[mnn_ka_dists.index]
    vals = np.ravel(lm.predict(x.values.reshape(-1, 1)))
    mnn_scaling_factors = pd.Series(vals, index=mnn_ka_dists.index)

    return mnn_scaling_factors


def _mnn_affinity(mnn, mnn_scaling_factors, offset1, offset2, device):
    # Function to convert mnn matrix to affinicty matrix

    # Construct adjacency matrix
    N = len(mnn_scaling_factors)
    x, y, z = find(mnn)
    x = x + offset1
    y = y + offset2
    adj = csr_matrix((z, (x, y)), shape=[N, N])

    # Affinity matrix
    return _convert_to_affinity(adj, mnn_scaling_factors, device, False)
