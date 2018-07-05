#!/usr/bin/env python


# Cluster_Ensembles/src/Cluster_Ensembles/Cluster_Ensembles.py;

# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com, ggiecold@jimmy.harvard.edu


"""Cluster_Ensembles is a package for combining multiple partitions 
into a consolidated clustering.
The combinatorial optimization problem of obtaining such a consensus clustering
is reformulated in terms of approximation algorithms for 
graph or hyper-graph partitioning.

References
----------
* Giecold, G., Marco, E., Trippa, L. and Yuan, G.-C.,
"Robust Lineage Reconstruction from High-Dimensional Single-Cell Data". 
ArXiv preprint [q-bio.QM, stat.AP, stat.CO, stat.ML]: http://arxiv.org/abs/1601.02748

* Strehl, A. and Ghosh, J., "Cluster Ensembles - A Knowledge Reuse Framework
for Combining Multiple Partitions".
In: Journal of Machine Learning Research, 3, pp. 583-617. 2002

* Kernighan, B. W. and Lin, S., "An Efficient Heuristic Procedure 
for Partitioning Graphs". 
In: The Bell System Technical Journal, 49, 2, pp. 291-307. 1970

* Karypis, G. and Kumar, V., "A Fast and High Quality Multilevel Scheme 
for Partitioning Irregular Graphs"
In: SIAM Journal on Scientific Computing, 20, 1, pp. 359-392. 1998

* Karypis, G., Aggarwal, R., Kumar, V. and Shekhar, S., "Multilevel Hypergraph Partitioning: 
Applications in the VLSI Domain".
In: IEEE Transactions on Very Large Scale Integration (VLSI) Systems, 7, 1, pp. 69-79. 1999
"""




import functools
import gc
import numbers
import numpy as np
import operator
import pkg_resources
import psutil
import scipy.sparse
from sklearn.metrics import jaccard_similarity_score
from sklearn.metrics import normalized_mutual_info_score
import subprocess
import sys
import tables
import warnings
import six
from six.moves import range
from functools import reduce

np.seterr(invalid = 'ignore')
warnings.filterwarnings('ignore', category = DeprecationWarning)


__all__ = ['cluster_ensembles', 'CSPA', 'HGPA', 'load_hypergraph_adjacency',
           'MCLA', 'overlap_matrix']


def memory():
    """Determine memory specifications of the machine.

    Returns
    -------
    mem_info : dictonary
        Holds the current values for the total, free and used memory of the system.
    """

    mem_info = dict()

    for k, v in psutil.virtual_memory()._asdict().items():
           mem_info[k] = int(v)
           
    return mem_info


def get_chunk_size(N, n):
    """Given a two-dimensional array with a dimension of size 'N', 
        determine the number of rows or columns that can fit into memory.

    Parameters
    ----------
    N : int
        The size of one of the dimensions of a two-dimensional array.  

    n : int
        The number of arrays of size 'N' times 'chunk_size' that can fit in memory.

    Returns
    -------
    chunk_size : int
        The size of the dimension orthogonal to the one of size 'N'. 
    """

    mem_free = memory()['free']
    if mem_free > 60000000:
        chunk_size = int(((mem_free - 10000000) * 1000) / (4 * n * N))
        return chunk_size
    elif mem_free > 40000000:
        chunk_size = int(((mem_free - 7000000) * 1000) / (4 * n * N))
        return chunk_size
    elif mem_free > 14000000:
        chunk_size = int(((mem_free - 2000000) * 1000) / (4 * n * N))
        return chunk_size
    elif mem_free > 8000000:
        chunk_size = int(((mem_free - 1400000) * 1000) / (4 * n * N))
        return chunk_size
    elif mem_free > 2000000:
        chunk_size = int(((mem_free - 900000) * 1000) / (4 * n * N))
        return chunk_size
    elif mem_free > 1000000:
        chunk_size = int(((mem_free - 400000) * 1000) / (4 * n * N))
        return chunk_size
    else:
        print("\nERROR: Cluster_Ensembles: get_chunk_size: "
              "this machine does not have enough free memory resources "
              "to perform ensemble clustering.\n")
        sys.exit(1)


def get_compression_filter(byte_counts):
    """Determine whether or not to use a compression on the array stored in
        a hierarchical data format, and which compression library to use to that purpose.
        Compression reduces the HDF5 file size and also helps improving I/O efficiency
        for large datasets.
    
    Parameters
    ----------
    byte_counts : int
    
    Returns
    -------
    FILTERS : instance of the tables.Filters class
    """

    assert isinstance(byte_counts, numbers.Integral) and byte_counts > 0
    
    if 2 * byte_counts > 1000 * memory()['free']:
        try:
            FILTERS = tables.filters(complevel = 5, complib = 'blosc', 
                                     shuffle = True, least_significant_digit = 6)
        except tables.FiltersWarning:
            FILTERS = tables.filters(complevel = 5, complib = 'lzo', 
                                     shuffle = True, least_significant_digit = 6)   
    else:
        FILTERS = None

    return FILTERS


def build_hypergraph_adjacency(cluster_runs):
    """Return the adjacency matrix to a hypergraph, in sparse matrix representation.
    
    Parameters
    ----------
    cluster_runs : array of shape (n_partitions, n_samples)
    
    Returns
    -------
    hypergraph_adjacency : compressed sparse row matrix
        Represents the hypergraph associated with an ensemble of partitions,
        each partition corresponding to a row of the array 'cluster_runs'
        provided at input.
    """

    N_runs = cluster_runs.shape[0]

    hypergraph_adjacency = create_membership_matrix(cluster_runs[0])
    for i in range(1, N_runs):
        hypergraph_adjacency = scipy.sparse.vstack([hypergraph_adjacency,
                                                   create_membership_matrix(cluster_runs[i])], 
                                                   format = 'csr')

    return hypergraph_adjacency


def store_hypergraph_adjacency(hypergraph_adjacency, hdf5_file_name):
    """Write an hypergraph adjacency to disk to disk in an HDF5 data structure.
    
    Parameters
    ----------
    hypergraph_adjacency : compressed sparse row matrix
    
    hdf5_file_name : file handle or string
    """
   
    assert(hypergraph_adjacency.__class__ == scipy.sparse.csr.csr_matrix)
    
    byte_counts = hypergraph_adjacency.data.nbytes + hypergraph_adjacency.indices.nbytes + hypergraph_adjacency.indptr.nbytes
    FILTERS = get_compression_filter(byte_counts)

    with tables.open_file(hdf5_file_name, 'r+') as fileh:
        for par in ('data', 'indices', 'indptr', 'shape'):
            try:
                n = getattr(fileh.root.consensus_group, par)
                n._f_remove()
            except AttributeError:
                pass

            array = np.array(getattr(hypergraph_adjacency, par))

            atom = tables.Atom.from_dtype(array.dtype)
            ds = fileh.create_carray(fileh.root.consensus_group, par, atom, 
                                     array.shape, filters = FILTERS)

            ds[:] = array


def load_hypergraph_adjacency(hdf5_file_name):
    """
    
    Parameters
    ----------
    hdf5_file_name : file handle or string
    
    Returns
    -------
    hypergraph_adjacency : compressed sparse row matrix
    """

    with tables.open_file(hdf5_file_name, 'r+') as fileh:
        pars = []
        for par in ('data', 'indices', 'indptr', 'shape'):
            pars.append(getattr(fileh.root.consensus_group, par).read())

    hypergraph_adjacency = scipy.sparse.csr_matrix(tuple(pars[:3]), shape = pars[3])
    
    return hypergraph_adjacency


def cluster_ensembles(cluster_runs, hdf5_file_name = None, verbose = False, N_clusters_max = None):
    """Call up to three different functions for heuristic ensemble clustering
       (namely CSPA, HGPA and MCLA) then select as the definitive
       consensus clustering the one with the highest average mutual information score 
       between its vector of consensus labels and the vectors of labels associated to each
       partition from the ensemble.

    Parameters
    ----------
    cluster_runs : array of shape (n_partitions, n_samples)
        Each row of this matrix is such that the i-th entry corresponds to the
        cluster ID to which the i-th sample of the data-set has been classified
        by this particular clustering. Samples not selected for clustering
        in a given round are are tagged by an NaN.
        
    hdf5_file_name : file object or string, optional (default = None)
        The handle or name of an HDF5 file where any array needed
        for consensus_clustering and too large to fit into memory 
        is to be stored. Created if not specified at input.
        
    verbose : Boolean, optional (default = False)
        Specifies if messages concerning the status of the many functions 
        subsequently called 'cluster_ensembles' will be displayed
        on the standard output.

    N_clusters_max : int, optional
        The number of clusters in which to partition the samples into 
        a consensus clustering. This defaults to the highest number of clusters
        encountered in the sets of independent clusterings on subsamples 
        of the data-set (i.e. the maximum of the entries in "cluster_runs").

    Returns
    -------
    cluster_ensemble : array of shape (n_samples,)
        For the final ensemble clustering, this vector contains the 
        cluster IDs of each sample in the whole data-set.

    Reference
    ---------
    A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse Framework
    for Combining Multiple Partitions".
    In: Journal of Machine Learning Research, 3, pp. 583-617. 2002  
    """

    if hdf5_file_name is None:
        hdf5_file_name = './Cluster_Ensembles.h5'
    fileh = tables.open_file(hdf5_file_name, 'w')
    fileh.create_group(fileh.root, 'consensus_group')
    fileh.close()

    cluster_ensemble = []
    score = np.empty(0)

    if cluster_runs.shape[1] > 10000:
        consensus_functions = [HGPA, MCLA]
        function_names = ['HGPA', 'MCLA']
        print("\nINFO: Cluster_Ensembles: cluster_ensembles: "
              "due to a rather large number of cells in your data-set, "
              "using only 'HyperGraph Partitioning Algorithm' (HGPA) "
              "and 'Meta-CLustering Algorithm' (MCLA) "
              "as ensemble consensus functions.\n")
    else:
        consensus_functions = [CSPA, HGPA, MCLA]
        function_names = ['CSPA', 'HGPA', 'MCLA']

    hypergraph_adjacency = build_hypergraph_adjacency(cluster_runs)
    store_hypergraph_adjacency(hypergraph_adjacency, hdf5_file_name)

    for i in range(len(consensus_functions)):
        cluster_ensemble.append(consensus_functions[i](hdf5_file_name, cluster_runs, verbose, N_clusters_max))
        score = np.append(score, ceEvalMutual(cluster_runs, cluster_ensemble[i], verbose))
        print("\nINFO: Cluster_Ensembles: cluster_ensembles: "
              "{0} at {1}.".format(function_names[i], score[i]))
        print('*****')

    return cluster_ensemble[np.argmax(score)]


def ceEvalMutual(cluster_runs, cluster_ensemble = None, verbose = False):
    """Compute a weighted average of the mutual information with the known labels, 
        the weights being proportional to the fraction of known labels.

    Parameters
    ----------
    cluster_runs : array of shape (n_partitions, n_samples)
        Each row of this matrix is such that the i-th entry corresponds to the
        cluster ID to which the i-th sample of the data-set has been classified
        by this particular clustering. Samples not selected for clustering
        in a given round are are tagged by an NaN.

    cluster_ensemble : array of shape (n_samples,), optional (default = None)
        The identity of the cluster to which each sample of the whole data-set 
        belong to according to consensus clustering.
 
    verbose : Boolean, optional (default = False)
        Specifies if status messages will be displayed
        on the standard output.

    Returns
    -------
    unnamed variable : float
        The weighted average of the mutual information between
        the consensus clustering and the many runs from the ensemble
        of independent clusterings on subsamples of the data-set.
    """

    if cluster_ensemble is None:
        return 0.0

    if reduce(operator.mul, cluster_runs.shape, 1) == max(cluster_runs.shape):
        cluster_runs = cluster_runs.reshape(1, -1)

    weighted_average_mutual_information = 0

    N_labelled_indices = 0

    for i in range(cluster_runs.shape[0]):
        labelled_indices = np.where(np.isfinite(cluster_runs[i]))[0]
        N = labelled_indices.size

        x = np.reshape(checkcl(cluster_ensemble[labelled_indices], verbose), newshape = N)
        y = np.reshape(checkcl(np.rint(cluster_runs[i, labelled_indices]), verbose), newshape = N)

        q = normalized_mutual_info_score(x, y)

        weighted_average_mutual_information += q * N
        N_labelled_indices += N

    return float(weighted_average_mutual_information) / N_labelled_indices


def checkcl(cluster_run, verbose = False):
    """Ensure that a cluster labelling is in a valid format. 

    Parameters
    ----------
    cluster_run : array of shape (n_samples,)
        A vector of cluster IDs for each of the samples selected for a given
        round of clustering. The samples not selected are labelled with NaN.

    verbose : Boolean, optional (default = False)
        Specifies if status messages will be displayed
        on the standard output.

    Returns
    -------
    cluster_run : array of shape (n_samples,)
        The input vector is modified in place, such that invalid values are
        either rejected or altered. In particular, the labelling of cluster IDs
        starts at zero and increases by 1 without any gap left.
    """
    
    cluster_run = np.asanyarray(cluster_run)

    if cluster_run.size == 0:
        raise ValueError("\nERROR: Cluster_Ensembles: checkcl: "
                         "empty vector provided as input.\n")
    elif reduce(operator.mul, cluster_run.shape, 1) != max(cluster_run.shape):
        raise ValueError("\nERROR: Cluster_Ensembles: checkl: "
                         "problem in dimensions of the cluster label vector "
                         "under consideration.\n")
    elif np.where(np.isnan(cluster_run))[0].size != 0:
        raise ValueError("\nERROR: Cluster_Ensembles: checkl: vector of cluster "
                         "labellings provided as input contains at least one 'NaN'.\n")
    else:
        min_label = np.amin(cluster_run)
        if min_label < 0:
            if verbose:
                print("\nINFO: Cluster_Ensembles: checkcl: detected negative values "
                      "as cluster labellings.")

            cluster_run -= min_label

            if verbose:
                print("\nINFO: Cluster_Ensembles: checkcl: "
                      "offset to a minimum value of '0'.")

        x = one_to_max(cluster_run) 
        if np.amax(cluster_run) != np.amax(x):
            if verbose:
                print("\nINFO: Cluster_Ensembles: checkcl: the vector cluster "
                      "labellings provided is not a dense integer mapping.")

            cluster_run = x

            if verbose:
                print("INFO: Cluster_Ensembles: checkcl: brought modification "
                      "to this vector so that its labels range "
                      "from 0 to {0}, included.\n".format(np.amax(cluster_run)))

    return cluster_run


def one_to_max(array_in):
    """Alter a vector of cluster labels to a dense mapping. 
        Given that this function is herein always called after passing 
        a vector to the function checkcl, one_to_max relies on the assumption 
        that cluster_run does not contain any NaN entries.

    Parameters
    ----------
    array_in : a list or one-dimensional array
        The list of cluster IDs to be processed.
    
    Returns
    -------
    result : one-dimensional array
        A massaged version of the input vector of cluster identities.
    """
    
    x = np.asanyarray(array_in)
    N_in = x.size
    array_in = x.reshape(N_in)    

    sorted_array = np.sort(array_in)
    sorting_indices = np.argsort(array_in)

    last = np.nan
    current_index = -1
    for i in range(N_in):
        if last != sorted_array[i] or np.isnan(last):
            last = sorted_array[i]
            current_index += 1

        sorted_array[i] = current_index

    result = np.empty(N_in, dtype = int)
    result[sorting_indices] = sorted_array

    return result


def checks(similarities, verbose = False):
    """Check that a matrix is a proper similarity matrix and bring 
        appropriate changes if applicable.

    Parameters
    ----------
    similarities : array of shape (n_samples, n_samples)
        A matrix of pairwise similarities between (sub)-samples of the data-set. 

    verbose : Boolean, optional (default = False)
        Alerts of any issue with the similarities matrix provided
        and of any step possibly taken to remediate such problem.
    """
    
    if similarities.size == 0:
        raise ValueError("\nERROR: Cluster_Ensembles: checks: the similarities "
                         "matrix provided as input happens to be empty.\n")
    elif np.where(np.isnan(similarities))[0].size != 0:
        raise ValueError("\nERROR: Cluster_Ensembles: checks: input similarities "
                         "matrix contains at least one 'NaN'.\n")
    elif np.where(np.isinf(similarities))[0].size != 0:
        raise ValueError("\nERROR: Cluster_Ensembles: checks: at least one infinite entry "
                         "detected in input similarities matrix.\n")
    else:
        if np.where(np.logical_not(np.isreal(similarities)))[0].size != 0:
            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: complex entries found "
                      "in the similarities matrix.")

            similarities = similarities.real

            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: "
                      "truncated to their real components.")

        if similarities.shape[0] != similarities.shape[1]:
            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: non-square matrix provided.")

            N_square = min(similarities.shape)
            similarities = similarities[:N_square, :N_square]

            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: using largest square sub-matrix.")

        max_sim = np.amax(similarities)
        min_sim = np.amin(similarities)
        if max_sim > 1 or min_sim < 0:
            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: strictly negative "
                      "or bigger than unity entries spotted in input similarities matrix.")

            indices_too_big = np.where(similarities > 1) 
            indices_negative = np.where(similarities < 0)
            similarities[indices_too_big] = 1.0
            similarities[indices_negative] = 0.0

            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: done setting them to "
                      "the lower or upper accepted values.")     

        if not np.allclose(similarities, np.transpose(similarities)):
            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: non-symmetric input "
                      "similarities matrix.")

            similarities = np.divide(similarities + np.transpose(similarities), 2.0)

            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: now symmetrized.")

        if not np.allclose(np.diag(similarities), np.ones(similarities.shape[0])):
            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: the self-similarities "
                      "provided as input are not all of unit value.")

            similarities[np.diag_indices(similarities.shape[0])] = 1

            if verbose:
                print("\nINFO: Cluster_Ensembles: checks: issue corrected.")


def CSPA(hdf5_file_name, cluster_runs, verbose = False, N_clusters_max = None):
    """Cluster-based Similarity Partitioning Algorithm for a consensus function.
    
    Parameters
    ----------
    hdf5_file_name : file handle or string
    
    cluster_runs : array of shape (n_partitions, n_samples)
    
    verbose : bool, optional (default = False)
    
    N_clusters_max : int, optional (default = None)
    
    Returns
    -------
    A vector specifying the cluster label to which each sample has been assigned
    by the CSPA heuristics for consensus clustering.

    Reference
    ---------
    A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse Framework
    for Combining Multiple Partitions".
    In: Journal of Machine Learning Research, 3, pp. 583-617. 2002
    """

    print('*****')
    print("INFO: Cluster_Ensembles: CSPA: consensus clustering using CSPA.")

    if N_clusters_max == None:
        N_clusters_max = int(np.nanmax(cluster_runs)) + 1

    N_runs = cluster_runs.shape[0]
    N_samples = cluster_runs.shape[1]
    if N_samples > 20000:
        raise ValueError("\nERROR: Cluster_Ensembles: CSPA: cannot efficiently "
                         "deal with too large a number of cells.")

    hypergraph_adjacency = load_hypergraph_adjacency(hdf5_file_name)

    s = scipy.sparse.csr_matrix.dot(hypergraph_adjacency.transpose().tocsr(), hypergraph_adjacency)
    s = np.squeeze(np.asarray(s.todense()))
    
    del hypergraph_adjacency
    gc.collect()

    checks(np.divide(s, float(N_runs)), verbose)

    e_sum_before = s.sum()
    sum_after = 100000000.0  
    scale_factor = sum_after / float(e_sum_before)

    with tables.open_file(hdf5_file_name, 'r+') as fileh:
        atom = tables.Float32Atom()
        FILTERS = get_compression_filter(4 * (N_samples ** 2))

        S = fileh.create_carray(fileh.root.consensus_group, 'similarities_CSPA', atom,
                               (N_samples, N_samples), "Matrix of similarities arising "
                               "in Cluster-based Similarity Partitioning", 
                               filters = FILTERS)

        expr = tables.Expr("s * scale_factor")
        expr.set_output(S)
        expr.eval()

        chunks_size = get_chunk_size(N_samples, 3)
        for i in range(0, N_samples, chunks_size):
            tmp = S[i:min(i+chunks_size, N_samples)]
            S[i:min(i+chunks_size, N_samples)] = np.rint(tmp)

    return metis(hdf5_file_name, N_clusters_max)


def HGPA(hdf5_file_name, cluster_runs, verbose = False, N_clusters_max = None):
    """HyperGraph-Partitioning Algorithm for a consensus function.
    
    Parameters
    ----------
    hdf5_file_name : string or file handle
    
    cluster_runs: array of shape (n_partitions, n_samples)
    
    verbose : bool, optional (default = False)
    
    N_clusters_max : int, optional (default = None)
    
    Returns
    -------
    A vector specifying the cluster label to which each sample has been assigned
    by the HGPA approximation algorithm for consensus clustering.

    Reference
    ---------
    A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse Framework
    for Combining Multiple Partitions".
    In: Journal of Machine Learning Research, 3, pp. 583-617. 2002
    """
    
    print('\n*****')
    print("INFO: Cluster_Ensembles: HGPA: consensus clustering using HGPA.")

    if N_clusters_max == None:
        N_clusters_max = int(np.nanmax(cluster_runs)) + 1

    return hmetis(hdf5_file_name, N_clusters_max)


def MCLA(hdf5_file_name, cluster_runs, verbose = False, N_clusters_max = None):
    """Meta-CLustering Algorithm for a consensus function.
    
    Parameters
    ----------
    hdf5_file_name : file handle or string
    
    cluster_runs : array of shape (n_partitions, n_samples)
    
    verbose : bool, optional (default = False)
    
    N_clusters_max : int, optional (default = None)
    
    Returns
    -------
    A vector specifying the cluster label to which each sample has been assigned
    by the MCLA approximation algorithm for consensus clustering.

    Reference
    ---------
    A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse Framework
    for Combining Multiple Partitions".
    In: Journal of Machine Learning Research, 3, pp. 583-617. 2002
    """

    print('\n*****')
    print('INFO: Cluster_Ensembles: MCLA: consensus clustering using MCLA.')

    if N_clusters_max == None:
        N_clusters_max = int(np.nanmax(cluster_runs)) + 1

    N_runs = cluster_runs.shape[0]
    N_samples = cluster_runs.shape[1]

    print("INFO: Cluster_Ensembles: MCLA: preparing graph for meta-clustering.")

    hypergraph_adjacency = load_hypergraph_adjacency(hdf5_file_name)
    w = hypergraph_adjacency.sum(axis = 1)

    N_rows = hypergraph_adjacency.shape[0]

    print("INFO: Cluster_Ensembles: MCLA: done filling hypergraph adjacency matrix. "
          "Starting computation of Jaccard similarity matrix.")

    # Next, obtain a matrix of pairwise Jaccard similarity scores between the rows of the hypergraph adjacency matrix.
    with tables.open_file(hdf5_file_name, 'r+') as fileh:
        FILTERS = get_compression_filter(4 * (N_rows ** 2))
    
        similarities_MCLA = fileh.create_carray(fileh.root.consensus_group, 
                                   'similarities_MCLA', tables.Float32Atom(), 
                                   (N_rows, N_rows), "Matrix of pairwise Jaccard "
                                   "similarity scores", filters = FILTERS)

        scale_factor = 100.0

        print("INFO: Cluster_Ensembles: MCLA: "
              "starting computation of Jaccard similarity matrix.")

        squared_MCLA = hypergraph_adjacency.dot(hypergraph_adjacency.transpose())

        squared_sums = hypergraph_adjacency.sum(axis = 1)
        squared_sums = np.squeeze(np.asarray(squared_sums))

        chunks_size = get_chunk_size(N_rows, 7)
        for i in range(0, N_rows, chunks_size):
            n_dim = min(chunks_size, N_rows - i)

            temp = squared_MCLA[i:min(i+chunks_size, N_rows), :].todense()
            temp = np.squeeze(np.asarray(temp))

            x = squared_sums[i:min(i+chunks_size, N_rows)]
            x = x.reshape(-1, 1)
            x = np.dot(x, np.ones((1, squared_sums.size)))

            y = np.dot(np.ones((n_dim, 1)), squared_sums.reshape(1, -1))
        
            temp = np.divide(temp, x + y - temp)
            temp *= scale_factor

            Jaccard_matrix = np.rint(temp)
            similarities_MCLA[i:min(i+chunks_size, N_rows)] = Jaccard_matrix

            del Jaccard_matrix, temp, x, y
            gc.collect()
 
    # Done computing the matrix of pairwise Jaccard similarity scores.
    print("INFO: Cluster_Ensembles: MCLA: done computing the matrix of "
          "pairwise Jaccard similarity scores.")

    cluster_labels = cmetis(hdf5_file_name, N_clusters_max, w)
    cluster_labels = one_to_max(cluster_labels)
    # After 'cmetis' returns, we are done with clustering hyper-edges

    # We are now ready to start the procedure meant to collapse meta-clusters.
    N_consensus = np.amax(cluster_labels) + 1

    fileh = tables.open_file(hdf5_file_name, 'r+')

    FILTERS = get_compression_filter(4 * N_consensus * N_samples)
    
    clb_cum = fileh.create_carray(fileh.root.consensus_group, 'clb_cum', 
                                  tables.Float32Atom(), (N_consensus, N_samples), 
                                  'Matrix of mean memberships, forming meta-clusters', 
                                  filters = FILTERS)  
 
    chunks_size = get_chunk_size(N_samples, 7)
    for i in range(0, N_consensus, chunks_size):
        x = min(chunks_size, N_consensus - i)
        matched_clusters = np.where(cluster_labels == np.reshape(np.arange(i, min(i + chunks_size, N_consensus)), newshape = (x, 1)))
        M = np.zeros((x, N_samples))
        for j in range(x):
            coord = np.where(matched_clusters[0] == j)[0]
            M[j] = np.asarray(hypergraph_adjacency[matched_clusters[1][coord], :].mean(axis = 0))
        clb_cum[i:min(i+chunks_size, N_consensus)] = M
    
    # Done with collapsing the hyper-edges into a single meta-hyper-edge, 
    # for each of the (N_consensus - 1) meta-clusters.

    del hypergraph_adjacency
    gc.collect()

    # Each object will now be assigned to its most associated meta-cluster.
    chunks_size = get_chunk_size(N_consensus, 4)
    N_chunks, remainder = divmod(N_samples, chunks_size)
    if N_chunks == 0:
        null_columns = np.where(clb_cum[:].sum(axis = 0) == 0)[0]
    else:
        szumsz = np.zeros(0)
        for i in range(N_chunks):
            M = clb_cum[:, i*chunks_size:(i+1)*chunks_size]
            szumsz = np.append(szumsz, M.sum(axis = 0))
        if remainder != 0:
            M = clb_cum[:, N_chunks*chunks_size:N_samples]
            szumsz = np.append(szumsz, M.sum(axis = 0))
        null_columns = np.where(szumsz == 0)[0]

    if null_columns.size != 0:
        print("INFO: Cluster_Ensembles: MCLA: {} objects with all zero associations "
              "in 'clb_cum' matrix of meta-clusters.".format(null_columns.size))
        clb_cum[:, null_columns] = np.random.rand(N_consensus, null_columns.size)

    random_state = np.random.RandomState()

    tmp = fileh.create_carray(fileh.root.consensus_group, 'tmp', tables.Float32Atom(),
                              (N_consensus, N_samples), "Temporary matrix to help with "
                              "collapsing to meta-hyper-edges", filters = FILTERS)

    chunks_size = get_chunk_size(N_samples, 2)
    N_chunks, remainder = divmod(N_consensus, chunks_size)
    if N_chunks == 0:
        tmp[:] = random_state.rand(N_consensus, N_samples)
    else:
        for i in range(N_chunks):
            tmp[i*chunks_size:(i+1)*chunks_size] = random_state.rand(chunks_size, N_samples)
        if remainder !=0:
            tmp[N_chunks*chunks_size:N_consensus] = random_state.rand(remainder, N_samples)

    expr = tables.Expr("clb_cum + (tmp / 10000)")
    expr.set_output(clb_cum)
    expr.eval()

    expr = tables.Expr("abs(tmp)")
    expr.set_output(tmp)
    expr.eval()

    chunks_size = get_chunk_size(N_consensus, 2)
    N_chunks, remainder = divmod(N_samples, chunks_size)
    if N_chunks == 0:
        sum_diag = tmp[:].sum(axis = 0)
    else:
        sum_diag = np.empty(0)
        for i in range(N_chunks):
            M = tmp[:, i*chunks_size:(i+1)*chunks_size]
            sum_diag = np.append(sum_diag, M.sum(axis = 0))
        if remainder != 0:
            M = tmp[:, N_chunks*chunks_size:N_samples]
            sum_diag = np.append(sum_diag, M.sum(axis = 0))

    fileh.remove_node(fileh.root.consensus_group, "tmp") 
    # The corresponding disk space will be freed after a call to 'fileh.close()'.

    inv_sum_diag = np.reciprocal(sum_diag.astype(float))

    if N_chunks == 0:
        clb_cum *= inv_sum_diag
        max_entries = np.amax(clb_cum, axis = 0)
    else:
        max_entries = np.zeros(N_samples)
        for i in range(N_chunks):
            clb_cum[:, i*chunks_size:(i+1)*chunks_size] *= inv_sum_diag[i*chunks_size:(i+1)*chunks_size]
            max_entries[i*chunks_size:(i+1)*chunks_size] = np.amax(clb_cum[:, i*chunks_size:(i+1)*chunks_size], axis = 0)
        if remainder != 0:
            clb_cum[:, N_chunks*chunks_size:N_samples] *= inv_sum_diag[N_chunks*chunks_size:N_samples]
            max_entries[N_chunks*chunks_size:N_samples] = np.amax(clb_cum[:, N_chunks*chunks_size:N_samples], axis = 0)

    cluster_labels = np.zeros(N_samples, dtype = int)
    winner_probabilities = np.zeros(N_samples)
    
    chunks_size = get_chunk_size(N_samples, 2)
    for i in reversed(range(0, N_consensus, chunks_size)):
        ind = np.where(np.tile(max_entries, (min(chunks_size, N_consensus - i), 1)) == clb_cum[i:min(i+chunks_size, N_consensus)])
        cluster_labels[ind[1]] = i + ind[0]
        winner_probabilities[ind[1]] = clb_cum[(ind[0] + i, ind[1])]       

    # Done with competing for objects.

    cluster_labels = one_to_max(cluster_labels)

    print("INFO: Cluster_Ensembles: MCLA: delivering "
          "{} clusters.".format(np.unique(cluster_labels).size))
    print("INFO: Cluster_Ensembles: MCLA: average posterior "
          "probability is {}".format(np.mean(winner_probabilities)))
    if cluster_labels.size <= 7:
        print("INFO: Cluster_Ensembles: MCLA: the winning posterior probabilities are:")
        print(winner_probabilities)
        print("'INFO: Cluster_Ensembles: MCLA: the full posterior probabilities are:")
        print(clb_cum)

    fileh.remove_node(fileh.root.consensus_group, "clb_cum")
    fileh.close()

    return cluster_labels


def create_membership_matrix(cluster_run):
    """For a label vector represented by cluster_run, constructs the binary 
        membership indicator matrix. Such matrices, when concatenated, contribute 
        to the adjacency matrix for a hypergraph representation of an 
        ensemble of clusterings.
    
    Parameters
    ----------
    cluster_run : array of shape (n_partitions, n_samples)
    
    Returns
    -------
    An adjacnecy matrix in compressed sparse row form.
    """

    cluster_run = np.asanyarray(cluster_run)

    if reduce(operator.mul, cluster_run.shape, 1) != max(cluster_run.shape):
        raise ValueError("\nERROR: Cluster_Ensembles: create_membership_matrix: "
                         "problem in dimensions of the cluster label vector "
                         "under consideration.")
    else:
        cluster_run = cluster_run.reshape(cluster_run.size)

        cluster_ids = np.unique(np.compress(np.isfinite(cluster_run), cluster_run))
      
        indices = np.empty(0, dtype = np.int32)
        indptr = np.zeros(1, dtype = np.int32)

        for elt in cluster_ids:
            indices = np.append(indices, np.where(cluster_run == elt)[0])
            indptr = np.append(indptr, indices.size)

        data = np.ones(indices.size, dtype = int)

        return scipy.sparse.csr_matrix((data, indices, indptr), shape = (cluster_ids.size, cluster_run.size))


def metis(hdf5_file_name, N_clusters_max):
    """METIS algorithm by Karypis and Kumar. Partitions the induced similarity graph 
        passed by CSPA.

    Parameters
    ----------
    hdf5_file_name : string or file handle
    
    N_clusters_max : int
    
    Returns
    -------
    labels : array of shape (n_samples,)
        A vector of labels denoting the cluster to which each sample has been assigned
        as a result of the CSPA heuristics for consensus clustering.
    
    Reference
    ---------
    G. Karypis and V. Kumar, "A Fast and High Quality Multilevel Scheme for
    Partitioning Irregular Graphs"
    In: SIAM Journal on Scientific Computing, Vol. 20, No. 1, pp. 359-392, 1999.
    """

    file_name = wgraph(hdf5_file_name)
    labels = sgraph(N_clusters_max, file_name)
    subprocess.call(['rm', file_name])

    return labels


def hmetis(hdf5_file_name, N_clusters_max, w = None):
    """Gives cluster labels ranging from 1 to N_clusters_max for 
        hypergraph partitioning required for HGPA.

    Parameters
    ----------
    hdf5_file_name : file handle or string
    
    N_clusters_max : int
    
    w : array, optional (default = None)
    
    Returns
    -------
    labels : array of shape (n_samples,)
        A vector of labels denoting the cluster to which each sample has been assigned
        as a result of the HGPA approximation algorithm for consensus clustering.
    
    Reference
    ---------
    G. Karypis, R. Aggarwal, V. Kumar and S. Shekhar, "Multilevel hypergraph
    partitioning: applications in VLSI domain" 
    In: IEEE Transactions on Very Large Scale Integration (VLSI) Systems, 
    Vol. 7, No. 1, pp. 69-79, 1999.
    """

    if w is None:
        file_name = wgraph(hdf5_file_name, None, 2)
    else:
        file_name = wgraph(hdf5_file_name, w, 3)
    labels = sgraph(N_clusters_max, file_name)
    labels = one_to_max(labels)

    subprocess.call(['rm', file_name])

    return labels


def cmetis(hdf5_file_name, N_clusters_max, w = None):
    """Returns cluster labellings ranging from 1 to N_clusters_max 
        for hypergraph partitioning involved in MCLA.

    Parameters
    ----------
    hdf5_file_name : file handle or string
    
    N_clusters_max : int
    
    w : array, optiona (default = None)

    Returns
    -------
    labels : array of shape (n_samples,)
        A vector of labels denoting the cluster to which each sample has been assigned
        as a result of the MCLA approximation algorithm for consensus clustering.
        
    Reference
    ---------
    G. Karypis and V. Kumar, "A Fast and High Quality Multilevel Scheme for
    Partitioning Irregular Graphs"
    In: SIAM Journal on Scientific Computing, Vol. 20, No. 1, pp. 359-392, 1999.
    """
 
    file_name = wgraph(hdf5_file_name, w, 1)
    labels = sgraph(N_clusters_max, file_name)
    labels = one_to_max(labels)

    subprocess.call(['rm', file_name])

    return labels


def wgraph(hdf5_file_name, w = None, method = 0):
    """Write a graph file in a format apposite to later use by METIS or HMETIS.
    
    Parameters
    ----------
    hdf5_file_name : file handle or string
    
    w : list or array, optional (default = None)
    
    method : int, optional (default = 0)
    
    Returns
    -------
    file_name : string
    """

    print('\n#')

    if method == 0:
        fileh = tables.open_file(hdf5_file_name, 'r+')
        e_mat = fileh.root.consensus_group.similarities_CSPA
        file_name = 'wgraph_CSPA'
    elif method == 1:
        fileh = tables.open_file(hdf5_file_name, 'r+')
        e_mat = fileh.root.consensus_group.similarities_MCLA
        file_name = 'wgraph_MCLA'
    elif method in {2, 3}:
        hypergraph_adjacency = load_hypergraph_adjacency(hdf5_file_name)
        e_mat = hypergraph_adjacency.copy().transpose()
        file_name = 'wgraph_HGPA'
        fileh = tables.open_file(hdf5_file_name, 'r+')
    else:
        raise ValueError("\nERROR: Cluster_Ensembles: wgraph: "
                         "invalid code for choice of method; "
                         "choose either 0, 1, 2 or 3.")

    if w is None:
        w = []

    N_rows = e_mat.shape[0]
    N_cols = e_mat.shape[1]

    if method in {0, 1}:
        diag_ind = np.diag_indices(N_rows)
        e_mat[diag_ind] = 0

    if method == 1:
        scale_factor = 100.0
        w_sum_before = np.sum(w)
        w *= scale_factor
        w = np.rint(w)

    with open(file_name, 'w') as file:
        print("INFO: Cluster_Ensembles: wgraph: writing {}.".format(file_name))

        if method == 0:
            sz = float(np.sum(e_mat[:] > 0)) / 2
            if int(sz) == 0:
                return 'DO_NOT_PROCESS'
            else:
                file.write('{} {} 1\n'.format(N_rows, int(sz)))
        elif method == 1:
            chunks_size = get_chunk_size(N_cols, 2)
            N_chunks, remainder = divmod(N_rows, chunks_size)
            if N_chunks == 0:
                sz = float(np.sum(e_mat[:] > 0)) / 2
            else:
                sz = 0
                for i in range(N_chunks):
                    M = e_mat[i*chunks_size:(i+1)*chunks_size]
                    sz += float(np.sum(M > 0))
                if remainder != 0:
                    M = e_mat[N_chunks*chunks_size:N_rows]
                    sz += float(np.sum(M > 0))
                sz = float(sz) / 2 
            file.write('{} {} 11\n'.format(N_rows, int(sz)))
        else:
            file.write('{} {} 1\n'.format(N_cols, N_rows))
                    
        if method in {0, 1}:
            chunks_size = get_chunk_size(N_cols, 2)
            for i in range(0, N_rows, chunks_size):
                M = e_mat[i:min(i+chunks_size, N_rows)]

                for j in range(M.shape[0]):
                    edges = np.where(M[j] > 0)[0]
                    weights = M[j, edges]

                    if method == 0:
                        interlaced = np.zeros(2 * edges.size, dtype = int)
                        # METIS and hMETIS have vertices numbering starting from 1:
                        interlaced[::2] = edges + 1 
                        interlaced[1::2] = weights
                    else:
                        interlaced = np.zeros(1 + 2 * edges.size, dtype = int)
                        interlaced[0] = w[i + j]
                        # METIS and hMETIS have vertices numbering starting from 1:
                        interlaced[1::2] = edges + 1 
                        interlaced[2::2] = weights

                    for elt in interlaced:
                        file.write('{} '.format(int(elt)))
                    file.write('\n')  
        else:
            print("INFO: Cluster_Ensembles: wgraph: {N_rows} vertices and {N_cols} "
                  "non-zero hyper-edges.".format(**locals()))

            chunks_size = get_chunk_size(N_rows, 2)
            for i in range(0, N_cols, chunks_size):
                M = np.asarray(e_mat[:, i:min(i+chunks_size, N_cols)].todense())
                for j in range(M.shape[1]):
                    edges = np.where(M[:, j] > 0)[0]
                    if method == 2:
                        weight = np.array(M[:, j].sum(), dtype = int)
                    else:
                        weight = w[i + j]
                    # METIS and hMETIS require vertices numbering starting from 1:
                    interlaced = np.append(weight, edges + 1) 
               
                    for elt in interlaced:
                        file.write('{} '.format(int(elt)))
                    file.write('\n')
    
    if method in {0, 1}:
        fileh.remove_node(fileh.root.consensus_group, e_mat.name)

    fileh.close()

    print('#')

    return file_name


def sgraph(N_clusters_max, file_name):
    """Runs METIS or hMETIS and returns the labels found by those 
        (hyper-)graph partitioning algorithms.
        
    Parameters
    ----------
    N_clusters_max : int
    
    file_name : string
    
    Returns
    -------
    labels : array of shape (n_samples,)
        A vector of labels denoting the cluster to which each sample has been assigned
        as a result of any of three approximation algorithms for consensus clustering 
        (either of CSPA, HGPA or MCLA).
    """

    if file_name == 'DO_NOT_PROCESS':
        return []

    print('\n#')

    k = str(N_clusters_max)
    out_name = file_name + '.part.' + k
    if file_name == 'wgraph_HGPA':
        print("INFO: Cluster_Ensembles: sgraph: "
              "calling shmetis for hypergraph partitioning.")
        
        if sys.platform.startswith('linux'):
            shmetis_path = pkg_resources.resource_filename(__name__, 
                                         'Hypergraph_Partitioning/hmetis-1.5-linux/shmetis')
        elif sys.platform.startswith('darwin'):
            shmetis_path = pkg_resources.resource_filename(__name__, 
                                      'Hypergraph_Partitioning/hmetis-1.5-osx-i686/shmetis')
        else:
            print("ERROR: Cluster_Ensembles: sgraph:\n"
                  "your platform is not supported. Some code required for graph partition "
                  "is only available for Linux distributions and OS X.")
            sys.exit(1)
        
        args = "{0} ./".format(shmetis_path) + file_name + " " + k + " 15"
        subprocess.call(args, shell = True)
    elif file_name == 'wgraph_CSPA' or file_name == 'wgraph_MCLA':
        print("INFO: Cluster_Ensembles: sgraph: "
              "calling gpmetis for graph partitioning.")
        args = "gpmetis ./" + file_name + " " + k
        subprocess.call(args, shell = True)
    else:
        raise NameError("ERROR: Cluster_Ensembles: sgraph: {} is not an acceptable "
                        "file-name.".format(file_name))

    labels = np.empty(0, dtype = int)
    with open(out_name, 'r') as file:
        print("INFO: Cluster_Ensembles: sgraph: (hyper)-graph partitioning completed; "
              "loading {}".format(out_name))
        labels = np.loadtxt(out_name, dtype = int)
        labels = labels.reshape(labels.size)
    labels = one_to_max(labels)            

    subprocess.call(['rm', out_name])

    print('#')

    return labels


def overlap_matrix(hdf5_file_name, consensus_labels, cluster_runs):
    """Writes on disk (in an HDF5 file whose handle is provided as the first
       argument to this function) a stack of matrices, each describing
       for a particular run the overlap of cluster ID's that are matching 
       each of the cluster ID's stored in 'consensus_labels' 
       (the vector of labels obtained by ensemble clustering). 
       Returns also the adjacency matrix for consensus clustering 
       and a vector of mutual informations between each of the clusterings 
       from the ensemble and their consensus.
       
    Parameters
    ----------
    hdf5_file_name : file handle or string
    
    consensus_labels : array of shape (n_samples,)
    
    cluster_runs : array of shape (n_partitions, n_samples)
    
    Returns
    -------
    cluster_dims_list : 
    
    mutual_info_list :
    
    consensus_adjacency :
    """

    if reduce(operator.mul, cluster_runs.shape, 1) == max(cluster_runs.shape):
        cluster_runs = cluster_runs.reshape(1, -1)

    N_runs, N_samples = cluster_runs.shape
    N_consensus_labels = np.unique(consensus_labels).size

    indices_consensus_adjacency = np.empty(0, dtype = np.int32)
    indptr_consensus_adjacency = np.zeros(1, dtype = np.int64)

    for k in range(N_consensus_labels):
        indices_consensus_adjacency = np.append(indices_consensus_adjacency, np.where(consensus_labels == k)[0])
        indptr_consensus_adjacency = np.append(indptr_consensus_adjacency, indices_consensus_adjacency.size)

    data_consensus_adjacency = np.ones(indices_consensus_adjacency.size, dtype = int) 

    consensus_adjacency = scipy.sparse.csr_matrix((data_consensus_adjacency, indices_consensus_adjacency, indptr_consensus_adjacency), 
                                                  shape = (N_consensus_labels, N_samples))

    fileh = tables.open_file(hdf5_file_name, 'r+')
    
    FILTERS = get_compression_filter(4 * N_consensus_labels * N_runs)

    overlap_matrix = fileh.create_earray(fileh.root.consensus_group, 'overlap_matrix',
                                         tables.Float32Atom(), (0, N_consensus_labels), 
                                         "Matrix of overlaps between each run and "
                                         "the consensus labellings", filters = FILTERS,
                                         expectedrows = N_consensus_labels * N_runs)

    mutual_info_list = []
    cluster_dims_list =  [0]

    for i in range(N_runs):
        M = cluster_runs[i]

        mutual_info_list.append(ceEvalMutual(M, consensus_labels))

        finite_indices = np.where(np.isfinite(M))[0]
        positive_indices = np.where(M >= 0)[0]
        selected_indices = np.intersect1d(finite_indices, positive_indices, assume_unique = True)
        cluster_ids = np.unique(M[selected_indices])
        n_ids = cluster_ids.size

        cluster_dims_list.append(n_ids)

        unions = np.zeros((n_ids, N_consensus_labels), dtype = float)

        indices = np.empty(0, dtype = int)
        indptr = [0]

        c = 0
        for elt in cluster_ids:
            indices = np.append(indices, np.where(M == elt)[0])
            indptr.append(indices.size)

            for k in range(N_consensus_labels):
                x = indices_consensus_adjacency[indptr_consensus_adjacency[k]:indptr_consensus_adjacency[k+1]]
                unions[c, k] = np.union1d(indices, x).size 
 
            c += 1 

        data = np.ones(indices.size, dtype = int)
    
        I = scipy.sparse.csr_matrix((data, indices, indptr), shape = (n_ids, N_samples))

        intersections = I.dot(consensus_adjacency.transpose())
        intersections = np.squeeze(np.asarray(intersections.todense()))

        overlap_matrix.append(np.divide(intersections, unions))

    fileh.close()

    return cluster_dims_list, mutual_info_list, consensus_adjacency
    
