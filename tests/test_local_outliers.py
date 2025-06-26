import pytest
import numpy as np
import pandas as pd
import anndata as ad
from spotsweeper.local_outliers import local_outliers
from spotsweeper.local_outliers import robust_z_score

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

def test_robust_z_score():
    """ check if robust z-score is calculated correctly """
    x = np.array([10, 1, 2, 3, 4])
    z = robust_z_score(x)
    expected = 0.6745 * (10 - 2.5) / 1.0
    assert pytest.approx(expected) == z

def test_local_outliers_cols(synthetic_adata):
    """test if the function adds expected columns"""
    adata = local_outliers(synthetic_adata, metric="detected", direction="both", n_neighbors=2)
    assert "detected_z" in adata.obs.columns
    assert "detected_outliers" in adata.obs.columns
    assert adata.obs["detected_outliers"].dtype == bool # check if it's a boolean type

def test_outlier_dim(synthetic_adata):
    """test if z-score and outliers have the correct dimension"""
    adata = local_outliers(synthetic_adata, metric="detected", direction="both", n_neighbors=2)
    assert len(adata.obs["detected_z"]) == synthetic_adata.n_obs # should be equal to # spots
    assert len(adata.obs["detected_outliers"]) == synthetic_adata.n_obs

def test_extreme_value(synthetic_adata):
    """test if extreme value is detected as an outlier"""
    adata = local_outliers(synthetic_adata, metric="detected", direction="higher", n_neighbors=2)
    sample_a_outliers = adata.obs.loc[adata.obs["sample_id"] == "A", "detected_outliers"].values # extract boolean array
    assert sample_a_outliers[-1]  # check if the value 300 has outlier indicator = True

def test_log_transform(synthetic_adata):
    """test if log transform creates corresponding columns"""
    adata = local_outliers(synthetic_adata, metric="detected", log=True, direction="lower", n_neighbors=2)
    assert "detected_log" in adata.obs
    assert "detected_z" in adata.obs

def test_invalid_metric(synthetic_adata):
    """test if invalid metric raises keyerror"""
    with pytest.raises(KeyError):
        local_outliers(synthetic_adata, metric="abcde")

def test_invalid_coord(synthetic_adata):
    """test if invalid coordinate key raises error"""
    with pytest.raises(KeyError):
        local_outliers(synthetic_adata, metric="detected", coord_key="abcde")