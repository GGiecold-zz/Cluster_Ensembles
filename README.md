# Cluster_Ensembles
A package for combining multiple partitions into a consolidated clustering. The combinatorial optimization problem of obtaining such a consensus clustering is reformulated in terms of approximation algorithms for graph or hyper-graph partitioning.

Installation
------------

Cluster_Ensembles is written in Python and in C. You need Python 2.7, its Standard Library and the following packages:
* NumPy (version 1.9.0 or any ulterior version);
* SciPy
* scikit-learn
* setuptools
* PyTables

As yet another prelimiary to running Cluster_Ensembles, you should also follow the instructions below:
* open a terminal console;
* type in ```sudo dnf install glibc.i686``` (for a machine running on Fedora) or ```sudo apt-get install libc6.i386`` (on a Ubuntu platform)

This will install the GNU C library that is required to run a 32-bit executable binary tasked with hyper-graph partitioning. Skipping this step would result in a 'bad ELF interpreter' error message when subsequently trying to run the Cluster_Ensembles package.

Upon completition of the steps outlined above, install Cluster_Ensembles by sending a request to the Python Package Index (PyPI) as follows:
* open a terminal console;
* enter ```pip install Cluster_Ensembles```

Any missing third-party dependency should be automatically resolved. As part of the installation of this package, some code written in C that will later on be required by the Cluster_Ensembles package to determine a graph partition is automatically compiled under the hood and according to the specifications of your machine. You therefore need to ensure availability of ```CMake``` and ```GNU make``` on your operating system.

Usage
-----

```
>>> import numpy as np
>>> import Cluster_Ensembles as CE
>>> cluster_runs = np.random.randint(0, 50, (50, 15000))
>>> consensus_clustering_labels = CE.cluster_ensembles(cluster_runs, verbose = True, N_clusters_max = 50)
```

References
----------
* Giecold, G., Marco, E., Trippa, L. and Yuan, G.-C.,
"Robust Inference of Cell Lineages", to appear
* A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse Framework
for Combining Multiple Partitions".
In: Journal of Machine Learning Research, 3, pp. 583-617. 2002

IMPORTANT NOTICE
----------------

A more detailed README file and expanded docstrings will be posted soon.
