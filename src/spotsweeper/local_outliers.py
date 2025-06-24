"""
local_outliers.py - Parallel to localOutliers.R
This function detects local outliers in spatial transcriptomics data based on
standard quality control metrics, such as library size, unique genes, and
mitochondrial ratio. Local outliers are defined as spots with low/high
quality metrics compared to their surrounding neighbors, based on a modified
z-score statistic.
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import anndata as ad


def local_outliers(
    adata: ad.AnnData,
    metric: str = "detected",
    direction: str = "lower",
    n_neighbors: int = 36,
    sample_key: str = "sample_id",
    log: bool = True,
    cutoff: float = 3.0,
    workers: int = 1,
    coord_key: str = "spatial",
) -> ad.AnnData:
    """
    This function detects local outliers in spatial transcriptomics data.

    Args:
    adata: AnnData object. MUST contain spatial coordinates in adata.obsm[coord_key] \
    and corresponding QC metric in adata.obs[metric]
    metric: column in adata.obs to use for outlier detection
    direction: direction of outlier detection (higher, lower, or both)
    n_neighbors: number of nearest neighbors to use for outlier detection
    sample_key: column name in adata.obs to use for sample IDs
    log: boolean indicating whether to log1p transform the metrics (default to TRUE)
    cutoff: cutoff for outlier detection (default = 3)
    workers: number of workers for parallel processing (default = 1)
    coord_key: key in adata.obsm containing spatial coordinates.

    Returns:
    adata: the updated AnnData object with updated columns containing outputs (in adata.obs)
    """

    # ----- Validity checks -----
    # check if coordinate key is in adata.obsm
    if coord_key not in adata.obsm:
        raise KeyError(
            f"Spatial coordinates '{coord_key}' not found in `adata.obsm`."
        )
    # check if metric is included in adata.obs
    if metric not in adata.obs:
        raise KeyError(f"Metric '{metric}' not found in `adata.obs`.")
    # validate direction
    if direction not in {"lower", "higher", "both"}:
        raise ValueError("`direction` must be 'lower', 'higher', or 'both'.")
    # check if n_neighbors is a positive integer
    if not (isinstance(n_neighbors, int) and n_neighbors > 0):
        raise ValueError("`n_neighbors` must be a positive integer.")
    # check if cutoff is a numeric value
    if not np.isscalar(cutoff):
        raise ValueError("`cutoff` must be a numeric value.")

    # log transform specified metric
    metric_to_use = metric
    if log:
        log_col = f"{metric}_log"
        if log_col not in adata.obs:  # if exists already, not overwrite
            adata.obs[log_col] = np.log1p(adata.obs[metric])
        metric_to_use = log_col

    # initialize result columns in AnnData
    z_col = f"{metric}_z"
    outlier_col = f"{metric}_outliers"
    adata.obs[z_col] = 0.0
    adata.obs[outlier_col] = False

    # get list of unique samples
    unique_samples = adata.obs[sample_key].unique()

    # loop through each unique sample ID
    for sample in unique_samples:
        # subset data for current sample
        sub_idx = adata.obs[sample_key] == sample
        adata_sub = adata[sub_idx]
        # extract coordinates and QC metrics
        coords = adata_sub.obsm[coord_key]
        values = adata_sub.obs[metric_to_use].to_numpy()

        # find nearest neighbors (train on spatial coordinates)
        # set the algorithm to "auto" select between brute force and tree-based methods
        nn = NearestNeighbors(
            n_neighbors=n_neighbors, algorithm="auto", n_jobs=workers
        ).fit(coords)
        neigh_idx = nn.kneighbors(return_distance=False) # neighbor index: 2D array with rows being spots, cols being neighbors

        # compute modified z for each spot
        z_scores = np.zeros(adata_sub.n_obs, dtype=float) # initialize 0s for # spots
        for i in range(adata_sub.n_obs):
            idx = np.concatenate(([i], neigh_idx[i])) # include self as the first element
            neighborhood = values[idx] # extract QC metrics
            z_scores[i] = robust_z_score(neighborhood)

        # handle non-finite z-scores and replace with 0
        z_scores[~np.isfinite(z_scores)] = 0.0

        # determine outliers based on cutoff
        if direction == "higher":
            outliers = z_scores > cutoff
        elif direction == "lower":
            outliers = z_scores < -cutoff # type: ignore
        else:  # both
            outliers = (z_scores > cutoff) | (z_scores < -cutoff) # type: ignore

        # write back to the original AnnData
        adata.obs.loc[sub_idx, z_col] = z_scores
        adata.obs.loc[sub_idx, outlier_col] = outliers

    return adata


def robust_z_score(x: np.ndarray) -> float:
    """
    Compute the modified/robust z-score of the first element in a 1D numpy array
    relative to other values in the array.

    Formula: 0.6745 * (xi - median(x)) / MAD(x),
    where MAD is the median absolute deviation.

    Args:
    x: 1D numpy array that includes metrics to be compared

    Returns:
    the robust z-score of the first element with respect to other elements.
    """
    if x.ndim != 1:
        raise ValueError("Input x must be a 1D array.") # validity check
    
    if len(x) <= 1:
        return 0.0 # no neighbors to compute the z-score

    # compute median and MAD of neighbors
    neighbors = x[1:]
    med = np.median(neighbors) 
    mad = np.median(np.abs(neighbors - med))
    
    if mad == 0 or not np.isfinite(mad): # handle bad values
        return 0.0
    
    return 0.6745 * (x[0] - med) / mad