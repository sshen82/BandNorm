# BandNorm

Ye Zheng\*, Siqi Shen\* and Sündüz Keleş. Normalization and De-noising of single-cell Hi-C Data with BandNorm and 3DVI. bioRxiv (2021). * contribute equally.

`BandNorm` is a Normalization method that removes depth effect and aligns distance effect efficiently. 
Users can input a path of multiple single cells, or a data.frame containing all cells and receive a normalized data.frame file.

<img src="https://github.com/sshen82/BandNorm/blob/main/figures/bandnorm_intro.png" alt="BandNorm" width="700px">

There are four functions, `download_schic`, `bandnorm`, `create_embedding`, and `plot_embedding`.
`download_schic` downloads one of the currently available single-cell Hi-C data cleaned by us,
`bandnorm` takes in sparse matrices and normalize them using BandNorm method, 
`create_embedding` summarizes the data into a PCA embedding in preparation for clustering and lower dimension embedding,
and `plot_embedding` calculates UMAP or tSNE embedding from the PCA obtained from `create_embedding`, and plots the resulting embedding.
(See vignette (`browseVignettes("BandNorm")) for more detail.)

## Installation

If necessary, install the dependencies:

```
install.packages(c('ggplot2', 'dplyr', 'data.table', 'Rtsne', 'umap'))

if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")

library(devtools)
install_github("immunogenomics/harmony")
```

`BandNorm` can be installed from Github:

```
devtools::install_github('sshen82/BandNorm', build_vignettes = TRUE)
library(BandNorm)
```

## Input

There are two possible ways to input the data:

1. A path containing all cells in the form of [chr1, binA, chr2, binB, count],
2. An R data.frame in the form of [chr1, binA, binB, count, diag, cell_name]. The column names should be c("chrom", "binA", "binB", "count", "diag", "cell")

## Download data

You can also use `download_schic` function to download real data. The available data includes Kim2020, Li2019, Ramani2017, and Lee2019. There are also summary files available for them, and the summary includes batch, cell-type, depth and sparsity information for them.

## Usage

```
data("hic_df")
data("batch")
data("cell_type")
bandnorm_result = bandnorm(hic_df = hic_df, save = FALSE)
embedding = create_embedding(hic_df = bandnorm_result, do_harmony = FALSE, batch = batch)
plot_embedding(embedding, "UMAP", cell_type = cell_type)
plot_embedding(embedding, "UMAP", cell_type = batch)
```
