"""
plot_QC.py - Parallel to plotQCmetrics.R
Plot QC metrics for a single sample in AnnData object
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import anndata as ad
from typing import Optional, Sequence # for type hints

def plot_qc_metrics(
    adata: ad.AnnData,
    sample_id: str = "sample_id",
    sample: Optional[str] = None,
    metric: str = "detected",
    outliers: Optional[str] = None,
    point_size: float = 2,
    colors: Sequence[str] = ("white", "black"),
    stroke: float = 1.0,
    coord_key: str = "spatial"
):
    """
    This function generates a plot for specified sample within AnnData object,
    highlighting outliers based on a specified metric. The plot visualizes the
    metric of interest and indicates outliers with a distinct color.

    Args:
    - adata: AnnData object. MUST contain spatial coordinates in adata.obsm[coord_key] \
    and corresponding QC metric in adata.obs[metric].
    - sample_id: a String that identifies column in adata.obs that contains sample IDs. Default to "sample_id".
    - sample: a String that identifies which sample to plot, default to None specified.
    - metric: a String that identifies the metric to be visualized. Must be a column name in adata.obs. \
    Default to "detected".
    - outliers: a String that specifies the column name in adata.obs that indicates \
    whether a data point is considered an outlier. Default to None.
    - point_size: a float value that specifies the size of points in the plot. Default to 2.
    - colors: a Sequence/list specifies colors to be used for gradient scale. \
    If length is 2, gradient will be single color gradient. Default to white,black
    - stroke: a float value that specifies border thickness for outlier points. Default to 1.
    - coord_key: key in adata.obsm containing spatial coordinates. Default to "spatial".

    Returns:
    plt: plot object created by matplotlib to visualize the specified metric and outliers. \
    The plot is not explicitly printed by the function and should be printed by the caller.
    """
    # subset adata to the specified sample
    if sample is None: # if no sample ID provided, default to the first unique ID
        sample = adata.obs[sample_id].unique()[0]
    mask = np.array(adata.obs[sample_id] == sample)
    adata_sub = adata[mask].copy() # copy to avoid issues

    # extract relevant data to build the plot
    coords = np.array(adata_sub.obsm[coord_key])
    df = pd.DataFrame(data=coords, index=adata_sub.obs_names, columns=["x","y"])
    df[metric] = adata_sub.obs[metric]

    # add outliers if they are present
    if outliers is not None:
        df["outlier"] = adata_sub.obs[outliers]
    else:
        df["outlier"] = False

    # build custom color scale
    from matplotlib.colors import LinearSegmentedColormap # to handle multiple colors
    if len(colors) == 2:
        cmap = sns.color_palette(palette=colors, as_cmap=True)
    elif len(colors) > 2:
        cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    else:
        raise ValueError("Color gradient must have at least 2 elements")

    # build plot
    plt.figure(figsize=(6, 6))
    scatter = plt.scatter(
        df["x"], df["y"], s = point_size**2, cmap = cmap,
        edgecolor=["red" if i else "none" for i in df["outlier"]], # red color for outliers
        linewidths=stroke
    )

    plt.title(f"Sample: {sample}")
    plt.axis("equal")
    plt.gca().invert_yaxis() # to match tissue orientation
    plt.colorbar(scatter, label=metric)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    return plt