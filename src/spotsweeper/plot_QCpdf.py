"""
plot_QCpdf.py - Parallel to plotQCpdf.R
Generate a PDF file containing plots for each sample in AnnData object,
highlighting outliers based on specified metrics. Each plot visualizes outlier
metrics for a single sample, allowing for easy comparison and analysis across samples.
"""
from matplotlib.backends.backend_pdf import PdfPages
from anndata import AnnData
from typing import Sequence
from plot_QC import plot_qc_metrics

def plot_qc_pdf(
    adata: AnnData,
    sample_id: str = "sample_id",
    metric: str = "detected",
    outliers: str = "local_outliers",
    colors: Sequence[str] = ("white", "black"),
    stroke: float = 1.0,
    point_size: float = 2.0,
    width: float = 5,
    height: float = 5,
    fname: str = "qc_plots.pdf"
):
    """
    Generate and save QC plots for each sample in AnnData object to PDF file.

    Args:
    - adata: AnnData object
    - sample_id: adata.obs column containing sample IDs. default to "sample_id"
    - metric: adata.obs column containing QC metric to plot. default to "detected"
    - outliers: obs column indicating outlier status. default to "local_outliers"
    - colors: color gradient for the metric. default to white, black
    - stroke: border thickness for outliers. default to 1.0
    - point_size: size of points in scatter plot. default to 2.0
    - width: width of plot in inches
    - height: height of plot in inches
    - fname: path and filename of the output PDF
    """
    sample_list = adata.obs[sample_id].unique()

    with PdfPages(fname) as pdf:
        for sample in sample_list:
            plt = plot_qc_metrics(
                adata=adata,
                sample_id=sample_id,
                sample=sample,
                metric=metric,
                outliers=outliers,
                point_size=point_size,
                colors=colors,
                stroke=stroke,   
            )
            fig = plt.gcf()
            fig.set_size_inches(width, height)
            pdf.savefig(fig)
            plt.close(fig)