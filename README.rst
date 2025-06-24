.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/spotsweeper_py.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/spotsweeper_py
    .. image:: https://readthedocs.org/projects/spotsweeper_py/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://spotsweeper_py.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/spotsweeper_py/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/spotsweeper_py
    .. image:: https://img.shields.io/pypi/v/spotsweeper_py.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/spotsweeper_py/
    .. image:: https://img.shields.io/conda/vn/conda-forge/spotsweeper_py.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/spotsweeper_py
    .. image:: https://pepy.tech/badge/spotsweeper_py/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/spotsweeper_py
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/spotsweeper_py

    .. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
        :alt: Project generated with PyScaffold
        :target: https://pyscaffold.org/


==============
SpotSweeper
==============


    Spatially-aware quality control for spatial transcriptomics


SpotSweeper is a package developed for spatially-aware quality control (QC) methods for the detection, visualization, and removal of both local outliers and regional artifacts in spot-based spatial transcriptomics data, such as 10x Genomics Visium, using standard QC metrics.

Features
--------
- Detect local outliers using robust z-score
- (to be added) Plotting functions of local outliers

Installation
------------

Install from PyPI:

.. code-block:: bash

    pip install spotsweeper

Usage
-----

.. code-block:: python

    import spotsweeper.local_outliers as lo 
    lo.local_outliers(adata, metric = "total_counts", sample_key = "region")

Project Status
--------------

This package is in early development. Use with caution; interfaces may change. 
Tests haven't been implemented and will be added in a future version.

.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.6. For details and usage
information on PyScaffold see https://pyscaffold.org/.
