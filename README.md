# Cluster_Ensembles
A package for combining multiple partitions into a consolidated clustering. The combinatorial optimization problem of obtaining such a consensus clustering is reformulated in terms of approximation algorithms for graph or hyper-graph partitioning.

Installation
------------

Cluster_Ensembles is written in Python and in C. You need Python 2.7, its Standard Library and the following packages:
* numexpr (version 2.4.0 or later)
* NumPy (version 1.9.0 or any ulterior version);
* SciPy
* scikit-learn
* setuptools
* PyTables

As yet another prelimiary to running Cluster_Ensembles, you should also follow the few more instructions below.

On CentOS, Fedora or some Red Hat Linux distribution:
* open a terminal console;
* type in: ```sudo dnf install glibc.i686```.

This will install the GNU C library that is required to run a 32-bit executable binary with a 64-bit Linux kernel. This executable is tasked with hyper-graph partitioning. Skipping this step would result in a ```bad ELF interpreter``` error message when subsequently trying to run the Cluster_Ensembles package.

On a Debian or Ubuntu platform, the following commands should yield the same outcome:
* open a terminal console;
* type in: ```sudo dpkg --add-architecture i386``` to add the i386 architecture;
* enter: ```sudo apt-get install libc6:i386```.

Upon completion of the steps outlined above, install Cluster_Ensembles by sending a request to the Python Package Index (PyPI) as follows:
* open a terminal console;
* enter ```pip install Cluster_Ensembles```.

Any missing third-party dependency should be automatically resolved. 
Please note that as part of the installation of this package, some code written in C that will later on be required by the Cluster_Ensembles package to determine a graph partition is automatically compiled under the hood and according to the specifications of your machine. 
You therefore need to ensure availability of ```CMake``` and ```GNU make``` on your operating system.

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
"Robust Inference of Cell Lineages", to submitted for publication
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

IMPORTANT NOTICE
----------------

A more detailed README file and expanded docstrings will be posted soon.
