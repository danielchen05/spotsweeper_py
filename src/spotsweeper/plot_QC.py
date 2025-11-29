"""
plot_QC.py - Parallel to plotQCmetrics.R
Plot QC metrics for a single sample in AnnData object
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
from typing import Optional, Sequence, Tuple  # added Tuple for figsize

def plot_qc_metrics(
    adata: ad.AnnData,
    sample_id: str = "sample_id",
    sample: Optional[str] = None,
    metric: str = "detected",
    outliers: Optional[str] = None,
    point_size: float = 2,
    colors: Sequence[str] = ("white", "black"),
    stroke: float = 1.0,
    coord_key: str = "spatial",
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (6, 6),
    legend: bool = False,   # optional legend
    ring_overlay: bool = True,  # 2-layer (Visium) vs 1-layer (Visium HD)
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
    - title: optional string for custom plot title. If None, defaults to f"Sample: {sample}".
    - figsize: tuple specifying figure size (width, height). Default to (6, 6).
    - legend: bool flag to add a legend distinguishing outliers vs non-outliers. Default to False.
    - ring_overlay: if True (default), use a 2-layer style where outliers are drawn as \
      red rings on top of a gradient background (recommended for standard Visium). \
      If False, use a 1-layer style where outliers have red edges on the same points \
      (recommended for dense Visium HD).

    Returns:
    plt: plot object created by matplotlib to visualize the specified metric and outliers. \
    The plot is not explicitly printed by the function and should be printed by the caller.
    """
    from matplotlib.colors import LinearSegmentedColormap  # to handle multiple colors
    from matplotlib.lines import Line2D

    # subset adata to the specified sample
    if sample is None:  # if no sample ID provided, default to the first unique ID
        sample = adata.obs[sample_id].unique()[0]
    mask = np.array(adata.obs[sample_id] == sample)
    adata_sub = adata[mask].copy()  # copy to avoid issues

    # extract relevant data to build the plot
    coords = np.array(adata_sub.obsm[coord_key])
    df = pd.DataFrame(data=coords, index=adata_sub.obs_names, columns=["x", "y"])
    df[metric] = adata_sub.obs[metric]

    # add outliers if they are present
    if outliers is not None:
        df["outlier"] = adata_sub.obs[outliers].astype(bool)
    else:
        df["outlier"] = False

    # build custom color scale
    if len(colors) >= 2:
        cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    else:
        raise ValueError("Color gradient must have at least 2 elements")

    # build plot
    plt.figure(figsize=figsize)  # custom control of figure size

    if ring_overlay:
        # 2-layer: base gradient + red rings (good for Visium)
        base = plt.scatter(
            df["x"],
            df["y"],
            c=df[metric],
            s=point_size**2,
            cmap=cmap,
            linewidths=0,
            rasterized=True,
        )

        # overlay: outliers as red rings
        highlighted = df[df["outlier"]]
        if outliers is not None and len(highlighted) > 0:
            plt.scatter(
                highlighted["x"],
                highlighted["y"],
                facecolors="none",
                edgecolors="red",
                linewidths=stroke * 1.5,
                s=(point_size * 1.4) ** 2,
            )
    else:
        # 1-layer: single scatter with red edges for outliers (better for Visium HD)
        base = plt.scatter(
            df["x"],
            df["y"],
            c=df[metric],
            s=point_size**2,
            cmap=cmap,
            edgecolors=["red" if i else "none" for i in df["outlier"]],
            linewidths=stroke,
        )

    plt.title(title if title is not None else f"Sample: {sample}")  # controlled title
    plt.axis("equal")
    plt.gca().invert_yaxis()  # to match tissue orientation
    plt.colorbar(base, label=metric)  # use base instead of scatter 
    plt.xlabel("x")
    plt.ylabel("y")

    # optional legend: outlier vs non-outlier
    if legend and outliers is not None:
        legend_elements = [
            Line2D([0], [0], marker="o", linestyle="None",
                markerfacecolor="white", markeredgecolor="black",
                markersize=6, label="Non-outlier"),
            Line2D([0], [0], marker="o", linestyle="None",
                markerfacecolor="white", markeredgecolor="red",
                markersize=6, label="Outlier"),
        ]
        plt.legend(handles=legend_elements, loc="upper right", frameon=True)

    plt.tight_layout()
    return plt