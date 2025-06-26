import pytest
import matplotlib
import numpy as np
import pandas as pd
import anndata as ad
matplotlib.use("Agg")
from spotsweeper.plot_QCpdf import plot_qc_pdf
from spotsweeper.local_outliers import local_outliers
from spotsweeper.plot_QC import plot_qc_metrics

@pytest.fixture # make sure synthetic data is passed in to every function
def synthetic_adata():
    """
    create synthetic data and return anndata object with coords, 2 samples, and QC metrics.
    """
    obs = pd.DataFrame({
        "sample_id": ["A"] * 5 + ["B"] * 5,
        "detected": [100, 105, 98, 102, 300, 50, 52, 48, 51, 49]
    })
    coords = np.array([
        [0, 0], [1, 0], [0, 1], [1, 1], [2, 2],  
        [0, 0], [1, 0], [0, 1], [1, 1], [2, 2]
    ])
    adata = ad.AnnData(obs=obs)
    adata.obsm["spatial"] = coords
    return adata

def test_pdf_exist(synthetic_adata, tmp_path):
    """tests if a valid PDF is created from a normal run"""
    adata = local_outliers(synthetic_adata,
                           metric="detected",
                           direction="higher",
                           n_neighbors=2)

    pdf_path = tmp_path / "qc.pdf"
    plot_qc_pdf(adata,
                sample_id="sample_id",
                metric="detected",
                outliers="detected_outliers",
                fname=str(pdf_path))
    assert pdf_path.exists()