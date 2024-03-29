# BandNorm

Ye Zheng\*, Siqi Shen\* and Sündüz Keleş. Normalization and De-noising of single-cell Hi-C Data with BandNorm and 3DVI. bioRxiv (2021). * contribute equally.

## What is BandNorm?

The advent of single-cell sequencing technologies in profiling 3D genome organization led to development of single-cell high-throughput chromatin conformation (scHi-C) assays. Data from these assays enhance our ability to study dynamic chromatin architecture and the impact of spatial genome interactions on cell regulation at an unprecedented resolution. At the individual cell resolution, heterogeneity driven by the stochastic nature of chromatin fiber, various nuclear processes, and unwanted variation due to sequencing depths and batch effects poses major analytical challenges for inferring single cell-level 3D genome organizations.

To explicitly capture chromatin conformation features and distinguish cells based on their 3D genome organizations, we develop a simple and fast band normalization approach, `BandNorm`, as well as a deep generative modeling framework, [3DVI](https://github.com/yezhengSTAT/3DVI), for more structured modeling of scHi-C data. `BandNorm` first removes genomic distance bias within a cell, and sequencing depth normalizes between cells. Consequently, `BandNorm` adds back a common band-dependent contact decay profile for the contact matrices across cells. The former step is achieved by dividing the interaction frequencies of each band within a cell with the cell’s band mean. The latter step is implemented by multiplying each scaled band by the average band mean across cells.

<img src="./figures/bandnorm_intro.png" alt="BandNorm" width="700px">

There are four functions, `download_schic`, `bandnorm`, `create_embedding`, and `plot_embedding`.
`download_schic` downloads one of the currently available single-cell Hi-C data cleaned by us,
`bandnorm` takes in sparse matrices and normalize them using BandNorm method, 
`create_embedding` summarizes the data into a PCA embedding in preparation for clustering and lower dimension embedding,
and `plot_embedding` calculates UMAP or tSNE embedding from the PCA obtained from `create_embedding`, and plots the resulting embedding.
(See vignette (`browseVignettes("BandNorm"))` or visit [BandNorm website](https://sshen82.github.io/BandNorm) for more detail.)

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

Also, if the user doesn't have default R path writing permission, it is possible to install the package to the personal library:

```
withr::with_libpaths("personalLibraryPath", devtools::install_github("sshen82/BandNorm"))
```

Below are the outputs of the code in the vignette.

<img src="./figures/cell_type.png" alt="celltype" width="700px">

<img src="./figures/batch.png" alt="batch" width="700px">

For more detail, please visit [BandNorm website](https://sshen82.github.io/BandNorm).
