# Cluster_Ensembles
A package for combining multiple partitions into a consolidated clustering. The combinatorial optimization problem of obtaining such a consensus clustering is reformulated in terms of approximation algorithms for graph or hyper-graph partitioning.

Installation
------------

Cluster_Ensembles is written in Python and in C. You need Python 2.7, its Standard Library and the following packages:
* numexpr (version 2.4.0 or later)
* NumPy (version 1.9.0 or any ulterior version)
* psutil
* SciPy
* scikit-learn
* setuptools
* PyTables

The ```pip install Cluster_Ensembles``` command mentioned below should automatically detect and, if applicable, install or update any of the afore-mentioned dependencies.

As yet another prelimiary to running Cluster_Ensembles, you should also follow the few more instructions below.

On CentOS, Fedora or some Red Hat Linux distribution:
* open a terminal console;
* type in: ```sudo dnf install glibc.i686```.

This will install the GNU C library that is required to run a 32-bit executable binary with a 64-bit Linux kernel. This executable is tasked with hyper-graph partitioning. Skipping this step would result in a ```bad ELF interpreter``` error message when subsequently trying to run the Cluster_Ensembles package.

On a Debian or Ubuntu platform, the following commands should yield the same outcome:
* open a terminal console;
* type in: ```sudo dpkg --add-architecture i386``` to add support for the i386 architecture;
* enter: ```sudo apt-get install libc6:i386```.

Upon completion of the steps outlined above, install Cluster_Ensembles by sending a request to the Python Package Index (PyPI) as follows:
* open a terminal console;
* enter ```pip install Cluster_Ensembles```.

Any missing third-party dependency should be automatically resolved. 
Please note that as part of the installation of this package, some code written in C that will later on be required by the Cluster_Ensembles package to determine a graph partition is automatically compiled under the hood and according to the specifications of your machine. 
You therefore need to ensure availability of ```CMake``` and ```GNU make``` on your operating system.

At any rate, to make sure that third-party dependecy software was automatically installed open a terminal window and type, for example, ```gpmetis -help```. If you don't get a help screen, please manually install [METIS - Serial Graph Partitioning and Fill-reducing Matrix Ordering] (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview). On Ubuntu run the command ```sudo apt-get install metis```.

Usage
-----

Say that you have an array of shape (M, N) where each row corresponds to a vector reporting the cluster label of each of the N samples comprising your dataset. It is possible that some of those samples have been left out of consideration from some of those M clusterings; in this case, the corresponding entry is tagged as NaN (```numpy.nan```). 

The few lines below illustrate how to submit consensus clustering analysis such an ```cluster_runs``` (M, N) array of cluster labels. 
A vector holding the consensus clustering identities for each of the N samples in your dataset, ```consensus_clustering_labels```, is returned.

Please note that those M vectors of clustering labels can correspond to partitions of the samples into distinct numbers of overall clusters. Cluster_Ensembles therefore offers the possibility of seeking a consensus clustering from the aggregation of a clustering of your dataset into, say, 10 groups, another clustering of a fraction of your samples into 5 clusters, yet another partition of your dataset into 20 clusters, etc. Those choices are entirely up to you. Pretty much all that is required for Cluster_Ensembles is an array of clustering vectors. 

```
>>> import numpy as np
>>> import Cluster_Ensembles as CE
>>> cluster_runs = np.random.randint(0, 50, (50, 15000))
>>> consensus_clustering_labels = CE.cluster_ensembles(cluster_runs, verbose = True, N_clusters_max = 50)
```

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
for Partitioning Irregular Graphs".
In: SIAM Journal on Scientific Computing, 20, 1, pp. 359-392. 1998
* Karypis, G., Aggarwal, R., Kumar, V. and Shekhar, S., "Multilevel Hypergraph Partitioning: 
Applications in the VLSI Domain".
In: IEEE Transactions on Very Large Scale Integration (VLSI) Systems, 7, 1, pp. 69-79. 1999
