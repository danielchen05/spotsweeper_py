# SpotSweeper

[![PyPI](https://img.shields.io/pypi/v/spotsweeper.svg)](https://pypi.org/project/spotsweeper/)
[![PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)

**Spatially-aware quality control for spatial transcriptomics**

SpotSweeper is a Python package developed for **spatially-aware quality control (QC)** of spot-based spatial transcriptomics data (e.g., 10x Genomics Visium and Visium HD).  
It enables detection and visualization of **local outliers** by comparing each spot’s QC metrics to its **local spatial neighborhood** using robust statistics.

---

## Manuscript

**Title:** *SpotSweeper-py: spatially-aware quality control metrics for spatial omics data in the Python ecosystem*  
**Authors:** Xingyi Chen, Michael Totty, Stephanie C. Hicks  
**Venue:** bioRxiv (2025)  
**DOI:** https://doi.org/10.64898/2025.12.06.692760  

If you use SpotSweeper-py, please cite the manuscript above.

---

## Features

- Detect **local outliers** using a modified / robust z-score
- Operate directly on `AnnData` objects
- Visualize QC metrics and outliers in spatial coordinates
- Export per-sample QC plots to multi-page PDF files
- Visualization styles suitable for both **Visium** and **Visium HD**

---

## Installation

Install from PyPI:

.. code-block:: bash

    pip install spotsweeper

Usage
-----

.. code-block:: python

    import spotsweeper.local_outliers as lo 
    import spotsweeper.plot_QC as plot_QC
    import spotsweeper.plot_QCpdf as pdf
    lo.local_outliers(adata, metric = "total_counts", sample_key = "region")
    plot_QC.plot_qc_metrics(adata,"region",metric = "total_counts", outliers="total_counts_outliers")
    pdf.plot_qc_pdf(adata,"region",metric = "total_counts", outliers="total_counts_outliers")

Project Status
--------------

This package is in early development. Use with caution; interfaces may change. 

.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.6. For details and usage
information on PyScaffold see https://pyscaffold.org/.
