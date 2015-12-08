#!/usr/bin/env python


# Cluster_Ensembles/src/Cluster_Ensembles/__init__.py;

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
Giecold, G., Marco, E., Trippa, L. and Yuan, G.-C., 
"Robust Inference of Cell Lineages", to appear

A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse Framework
for Combining Multiple Partitions".
In: Journal of Machine Learning Research, 3, pp. 583-617. 2002
"""


from . Cluster_Ensembles import *


__all__ = []
__all__ += Cluster_Ensembles.__all__

