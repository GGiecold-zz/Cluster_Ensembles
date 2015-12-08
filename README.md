# Cluster_Ensembles
A package for combining multiple partitions into a consolidated clustering. The combinatorial optimization problem of obtaining such a consensus clustering is reformulated in terms of approximation algorithms for graph or hyper-graph partitioning.

Installation
------------

Cluster_Ensembles involves some Python and C code. You need Python 2.7, its Standard Library and the following packages:
* NumPy (versioon 1.9.0 or any ulterior version);
* SciPy
* scikit-learn
* setuptools
* PyTables

As yet another prelimiary to running Cluster_Ensembles, you should also follow the instructions below:
* open a terminal console;
* type in ```sudo dnf install glibc.i686``` (for a machine running on Fedora) or ```sudo apt-get install libc6.i386`` (on a Ubuntu platform)

This will install the GNU C library required to run a 32-bit executable binary tasked with hyper-graph partitioning. Skipping this step will later account for an 'bad ELF interpreter' error message when trying to run the Cluster_Ensembles package.

Upon completition of the steps outlined above, install Cluster_Ensembles by sending a request to the Python Package Index (PyPI) as follows:
* open a terminal console;
* enter ```pip install Cluster_Ensembles```

Any missing third-party dependency should be automatically resolved. The installation of this package also automatically compile under the hood and according to the specifications of your machine some code written in C that will later on be required by the Cluster_Ensembles package so as to determine a graph partition.

References
----------
* Giecold, G., Marco, E., Trippa, L. and Yuan, G.-C.,
"Robust Inference of Cell Lineages", to appear
* A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse Framework
for Combining Multiple Partitions".
In: Journal of Machine Learning Research, 3, pp. 583-617. 2002

Please Note
-----------

A more detailed README file and expanded docstrings will be posted soon.
