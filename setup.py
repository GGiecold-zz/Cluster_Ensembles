#!/usr/bin/env python


# Cluster_Ensembles/setup.py;

# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com, ggiecold@jimmy.harvard.edu


"""Setup script for Cluster_Ensembles, a package for combining multiple
partitions into a consolidated clustering.
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


from codecs import open
from os import path
from sys import exit, version
from setuptools import setup
from setuptools.command.install import install
import subprocess


here = path.abspath(path.dirname(__file__))


class My_install(install):

    def run(self):
        try:    
            subprocess.call(['make config'], cwd = path.join(here, 'src', 
                            'Cluster_Ensembles', 'Hypergraph_Partitioning', 'metis-5.1.0'),
                            shell = True)
            subprocess.call(['make'], cwd = path.join(here, 'src',
                            'Cluster_Ensembles', 'Hypergraph_Partitioning', 'metis-5.1.0'),
                            shell = True)
        except Exception as e:
            print(e)
            print("ERROR: Cluster_Ensembles: setup:\n"
                  "error occurred while attempting to compile metis: "
                  "try running 'make' instead.")
            exit(1)
        else:
            install.run(self)


with open(path.join(here, 'README.md'), encoding = 'utf-8') as f:
    long_description = f.read()
    

setup(name = 'Cluster_Ensembles',
      version = '1.16',
      
      description = "A package for determining the consensus clustering from " 
                    "an ensemble of partitions",
      long_description = long_description,
                    
      url = 'https://github.com/GGiecold/Cluster_Ensembles',
      download_url = 'https://github.com/GGiecold/Cluster_Ensembles',
      
      author = 'Gregory Giecold',
      author_email = 'g.giecold@gmail.com',
      maintainer = 'Gregory Giecold',
      maintainer_email = 'ggiecold@jimmy.harvard.edu',
      
      license = 'MIT License',
      
      platforms = ('Any',),
      install_requires = ['numexpr>=2.4', 'numpy>=1.9.0', 'psutil', 'scipy', 'sklearn', 
                          'setuptools', 'tables'],
                          
      classifiers = ['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: End Users/Desktop',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',          
                   'License :: OSI Approved :: MIT License',
                   'Natural Language :: English',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: POSIX',
                   'Programming Language :: C',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Scientific/Engineering :: Visualization',
                   'Topic :: Scientific/Engineering :: Mathematics', ],
                   
      packages = ['Cluster_Ensembles'],
      package_dir = {'Cluster_Ensembles': 'src/Cluster_Ensembles'},
      
      include_package_data = True,
      package_data = {
          'Cluster_Ensembles': 
              ['Hypergraph_Partitioning/hmetis-1.5-linux/hmetis',
               'Hypergraph_Partitioning/hmetis-1.5-linux/khmetis',
               'Hypergraph_Partitioning/hmetis-1.5-linux/shmetis',
               'Hypergraph_Partitioning/hmetis-1.5-osx-i686/hmetis',
               'Hypergraph_Partitioning/hmetis-1.5-osx-i686/khmetis',
               'Hypergraph_Partitioning/hmetis-1.5-osx-i686/shmetis',
               'Hypergraph_Partitioning/metis-5.1.0/*.txt',
               'Hypergraph_Partitioning/metis-5.1.0/Changelog',
               'Hypergraph_Partitioning/metis-5.1.0/Makefile',
               'Hypergraph_Partitioning/metis-5.1.0/vsgen.bat',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/*.c',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/*.h',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/*.txt',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/GKlibSystem.cmake',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/Makefile',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/conf/check_thread_storage.c',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/test/*.c',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/test/*.txt',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/test/Makefile.in.old',
               'Hypergraph_Partitioning/metis-5.1.0/GKlib/test/Makefile.old',
               'Hypergraph_Partitioning/metis-5.1.0/graphs/4elt.graph',
               'Hypergraph_Partitioning/metis-5.1.0/graphs/copter2.graph',
               'Hypergraph_Partitioning/metis-5.1.0/graphs/mdual.graph',
               'Hypergraph_Partitioning/metis-5.1.0/graphs/metis.mesh',
               'Hypergraph_Partitioning/metis-5.1.0/graphs/README',
               'Hypergraph_Partitioning/metis-5.1.0/graphs/test.mgraph',
               'Hypergraph_Partitioning/metis-5.1.0/include/CMakeLists.txt',
               'Hypergraph_Partitioning/metis-5.1.0/include/metis.h',
               'Hypergraph_Partitioning/metis-5.1.0/libmetis/*.c',
               'Hypergraph_Partitioning/metis-5.1.0/libmetis/CMakeLists.txt',
               'Hypergraph_Partitioning/metis-5.1.0/libmetis/*.h',
               'Hypergraph_Partitioning/metis-5.1.0/programs/CMakeLists.txt',
               'Hypergraph_Partitioning/metis-5.1.0/programs/*.c',
               'Hypergraph_Partitioning/metis-5.1.0/programs/*.h'          
              ],
      },
      
      cmdclass = {'install': My_install},
                   
      keywords = "aggregation clustering consensus consensus-clustering CSPA "
                 "data-mining ensemble ensemble-clustering HGPA hyper-graph "
                 "machine-learning MCLA partition pattern-recognition "
                 "unsupervised-learning", 
)

