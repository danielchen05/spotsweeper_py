# SpotSweeper-py

[![PyPI-Server](https://img.shields.io/pypi/v/spotsweeper.svg)](https://pypi.org/project/spotsweeper/)
[![Project generated with PyScaffold](https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold)](https://pyscaffold.org/)

Spatially-aware quality control for spatial transcriptomics

SpotSweeper-py is a PyPI package developed for spatially-aware quality control (QC)
methods for the detection, visualization, and removal of local outliers in spot-based spatial transcriptomics data (e.g., 10x Genomics
Visium and Visium HD), using standard QC metrics.

---

## Manuscript

**Title:** *SpotSweeper-py: spatially-aware quality control metrics for spatial omics data in the Python ecosystem*  
**Authors:** Xingyi Chen, Michael Totty, Stephanie C. Hicks  
**Venue:** bioRxiv (2025)  
**DOI:** https://doi.org/10.64898/2025.12.06.692760

If you use SpotSweeper-py, please cite the manuscript above.

---

## Features

- Detect local outliers using a modified / robust z-score
- Operate directly on `AnnData` objects
- Visualize QC metrics and highlight local outliers in spatial context
- Export per-sample QC plots to multi-page PDF files
- Visualization styles suitable for both Visium and Visium HD

---

## Installation

Install from PyPI:

```bash
pip install spotsweeper
```

---

## Dependencies

SpotSweeper-py depends on the following core Python packages:

- `numpy`
- `pandas`
- `scikit-learn`
- `anndata`
- `matplotlib`

All required dependencies are installed automatically when installing
SpotSweeper-py from PyPI.

---

## Usage

SpotSweeper-py operates directly on `AnnData` objects. A typical workflow is:

- Detect local outliers for a QC metric (results are written to `adata.obs`)
- Visualize a QC metric for a single sample, optionally highlighting outliers
- Optionally export per-sample QC plots to a multi-page PDF

---

## Quickstart (local outliers + plot + optional PDF)

```python
import numpy as np
import spotsweeper.local_outliers as lo
import spotsweeper.plot_QC as plot_QC
import spotsweeper.plot_QCpdf as pdf

# Compute log total counts (common QC transform)
# Skip this step if the column already exists
adata.obs["log_total_counts"] = np.log1p(adata.obs["total_counts"])

# Detect local outliers using log total counts
lo.local_outliers(
    adata,
    metric="log_total_counts",
    direction="lower",
    n_neighbors=36,
    sample_key="region",
    log=False,
    cutoff=3.0,
    coord_key="spatial",
)

# Visualize local outliers for a single sample
plot_QC.plot_qc_metrics(
    adata,
    "region",
    metric="log_total_counts",
    outliers="log_total_counts_outliers",
    title="SpotSweeper QC",
    legend=True,
)

# (Optional) Save per-sample QC plots to a PDF
pdf.plot_qc_pdf(
    adata,
    "region",
    metric="log_total_counts",
    outliers="log_total_counts_outliers",
    fname="qc_plots.pdf",
)
```

---

## Common QC metrics

SpotSweeper-py can be applied to any numeric QC metric stored in `adata.obs`.
Common choices include:

- `total_counts` (library size / UMI count)
- `log_total_counts` (log-transformed library size)
- `n_genes_by_counts` (number of detected genes)
- `pct_counts_mt` (mitochondrial fraction)

If a raw metric (e.g., `total_counts`) is supplied and `log=True` (default),
SpotSweeper-py will internally apply `log1p` and store the transformed values as
`<metric>_log` before computing local z-scores.

If a precomputed metric (e.g., `log_total_counts`) is supplied, set `log=False`
to avoid double transformation.

---

## Choosing the outlier direction

Use the `direction` argument to control which tail(s) are flagged.
By default, `direction="lower"`.

- `direction="lower"` (default): flags unusually low metric values  (e.g., low counts or low numbers of genes)
- `direction="higher"`: flags unusually high metric values  (e.g., high mitochondrial fraction)
- `direction="both"`: flags both tails

---

## Example (high mitochondrial fraction)

```python
import spotsweeper.local_outliers as lo

lo.local_outliers(
    adata,
    metric="pct_counts_mt",
    direction="higher",
    sample_key="sample_id",
    cutoff=3.0,
)
```

---

## Plot styling for Visium vs Visium HD

The plotting function supports two visualization styles via `ring_overlay`:

- `ring_overlay=True` (default): two-layer plot with metric gradient and red rings
  for outliers; recommended for standard Visium data.
- `ring_overlay=False`: single-layer plot with red edges for outliers; recommended
  for dense Visium HD data.

---

## Example (dense data; single-layer style)

```python
import spotsweeper.plot_QC as plot_QC

plot_QC.plot_qc_metrics(
    adata,
    sample_id="sample_id",
    metric="detected",
    outliers="detected_outliers",
    ring_overlay=False,
    legend=True,
)
```

---

## Example notebooks

End-to-end example notebooks reproducing analyses and figures from the manuscript
are available in the companion analysis repository:

https://github.com/danielchen05/SpotSweeper_py_paper

---

## Requirements

SpotSweeper-py expects the input `AnnData` object to contain:

- Spatial coordinates stored in `adata.obsm["spatial"]`
  (or another key specified by `coord_key`)
- QC metrics stored as columns in `adata.obs`
  (e.g., `total_counts`, `detected`, `pct_counts_mt`)
- A sample identifier column in `adata.obs`
  (e.g., `sample_id` or `region`)

The package does not perform data loading or preprocessing and is agnostic to how
the `AnnData` object is constructed.

---

## Project Status

SpotSweeper-py is a PyPI software package accompanying a bioRxiv preprint.
The core methodology and functionality are stable and documented in the manuscript,
while the software interface may continue to evolve with additional features and
improvements.

---

## Contributing

Bug reports, feature requests, and GitHub issues or pull requests are welcome.

Please submit issues and pull requests via the GitHub repository:
https://github.com/danielchen05/spotsweeper_py

---

## Note

This project has been set up using PyScaffold 4.6.
For details and usage information on PyScaffold see https://pyscaffold.org/.