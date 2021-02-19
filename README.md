# BandNorm

`BandNorm` is a Normalization method that removes depth effect and aligns distance effect efficiently. 
Users can input a path of multiple single cells, or a data.frame containing all cells and receive a normalized data.frame file.

There are three functions, `bandnorm`, `create_embedding`, and `plot_embedding`.
`bandnorm` takes in sparse matrices and normalize them using BandNorm method, 
`create_embedding` summarizes the data into a PCA embedding in preparation for clustering and lower dimension embedding,
and `plot_embedding` calculates UMAP or tSNE embedding from the PCA obtained from `create_embedding`, and plots the resulting embedding.
(See vignette (`browseVignettes("BandNorm")) for more detail.)

## Installation

If necessary, install the dependencies:

```
install.packages(c('ggplot2', 'dplyr', 'data.table',
    'Rtsne',
    'umap'))

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
