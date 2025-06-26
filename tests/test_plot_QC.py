import pytest
import matplotlib
import numpy as np
import pandas as pd
import anndata as ad
matplotlib.use("Agg") # do not display
import matplotlib.pyplot as plt
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

def test_plot_axis(synthetic_adata):
    """
    test if plot function creates a figure, y-axis is inverted, and color-bar successfully added
    """
    plot_qc_metrics(
        synthetic_adata,
        sample_id="sample_id",
        sample="A",
        metric="detected",
        colors=("white", "black"),
        point_size=3,
        coord_key="spatial",
    )

    fig = plt.gcf()
    ax  = fig.axes[0]
    assert ax.yaxis_inverted()
    assert len(fig.axes) >= 2 # colorbar has its own axis
    plt.close(fig)

def test_color_list(synthetic_adata):
    """test if passing one color will raise error"""
    with pytest.raises(ValueError):
        plot_qc_metrics(
            synthetic_adata,
            colors=("blue",),  # only one colour supplied
            sample="A",
            metric="detected",
        )
