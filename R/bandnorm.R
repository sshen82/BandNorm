#' BandNorm
#'
#' This function allows you to calculate the BandNorm normalization.
#' @param path The path for all the cells in a directory. There can be sub-directories.
#' @param hic_df If you prepare the dataset as a format of "chrom", "binA", "binB", "count", "diag", "cell", you can input it into bandnorm directly.
#' @param save Whether to save each normalized cells. Default is TRUE. Note that if don't have large memory on your computer, and you need to use create_embedding function, it is highly recommended to save the cells because it helps lower the cost of memory in this function.
#' @param save_path Indicate the output path for normalized cells. Only need it when "save" parameter is TRUE. Default is NULL.
#' @export
#' @import data.table
#' @import dplyr
#' @examples
#' data("hic_df")
#' bandnorm_result = bandnorm(hic_df = hic_df, save = FALSE)
bandnorm = function(path = NULL, hic_df = NULL, save = TRUE, save_path = NULL) {
  # Get path and name for all the cells in this path.
  if (!is.null(path)){
    paths = list.files(path, recursive = TRUE, full.names = TRUE)
    names = list.files(path, recursive = TRUE)
    # The input format of the cell should be [chr1, bin1, chr2, bin2, count].
    load_cell = function(i) {
      return(fread(paths[i]) %>% rename(chrom = V1, binA = V2, binB = V4, count = V5) %>%
               mutate(diag = abs(binB - binA), cell = names[i]) %>% select(-V3))
    }
    hic_df = rbindlist(lapply(1:length(paths), load_cell))
  }
  # Calculate band depth and the mean of band depth for bandnorm.
  band_info <- hic_df %>% group_by(chrom, diag, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(chrom, diag) %>% summarise(depth = mean(band_depth))

  hic_df <- hic_df %>% left_join(alpha_j, by = c("chrom", "diag")) %>% left_join(band_info,
                                                                                 by = c("chrom", "diag", "cell")) %>% mutate(BandNorm = ifelse(diag == 0,
                                                                                 count, count/band_depth * depth)) %>% select(-c(band_depth, depth, count))
  if (save) {
    save_path = file.path(save_path, unique(names))
    for (i in 1:length(unique(names))) {
      write.table(hic_df[cell == unique(names)[i], c("chrom", "binA", "chrom",
                                                     "binB", "BandNorm")], file = save_path[i], col.names = FALSE, row.names = FALSE,
                  quote = FALSE, sep = "\t")
    }
  }
  return(hic_df)
}

#' create_embedding
#'
#' This function allows you to obtain the PCA embedding for tSNE or UMAP plotting.
#' @param path The path for all the normalized cells in a directory. There can be sub-directories. Using path means loading the bandnorm normalized data iteratively from the directory, so it is relatively slower than using hic_df. However, it won't eat up your memory too much, and the speed is acceptable.
#' @param hic_df After using bandnorm, if you keep the data frame, it is possible to use this instead of "path" as input of create_embedding. If using hic_df, the speed will be faster, but is costs more memory.
#' @param mean_thres The proportion of low mean interactions to be filtered. Ranges from 0 to 1, default is 0.
#' @param var_thres The proportion of low variance interactions to be filtered. Ranges from 0 to 1, default is 0.
#' @param dim_pca Dimension of PCA embedding to be outputted. Default is 50.
#' @param do_harmony Whether to use Harmony to remove the batch effect from the embedding. Default is FALSE
#' @param batch The batch information used for Harmony to remove the batch effect. Required if do_harmony is TRUE.
#' @export
#' @import data.table
#' @import dplyr
#' @import harmony
#' @examples
#' data("hic_df")
#' data("batch")
#' bandnorm_result = bandnorm(hic_df = hic_df, save = FALSE)
#' embedding = create_embedding(hic_df = bandnorm_result, do_harmony = TRUE, batch = batch)
create_embedding = function(path = NULL, hic_df = NULL, mean_thres = 0, var_thres = 0,
                            dim_pca = 50, do_harmony = FALSE, batch = NULL) {
  # Function to create embedding for cells, after combining all the bin-pairs from all chromosomes,
  # we do PCA first, and this function returns the PCA embedding.

  # There are two options:
  # If using hic_df, the speed will be faster, but is costs more memory, so you can use it if the memory
  # of the computer is sufficient;
  # For path, it means loading the bandnorm normalized data iteratively, so it is slower. However,
  # it won't eat up your memory too much, and the speed is also acceptable.
  if (is.null(path)) {
    setDT(hic_df)
    summarized_hic = hic_df[, .(agg_m = mean(BandNorm), agg_v = var(BandNorm)),
                            by = .(chrom, binA, binB)]
    summarized_hic[is.na(agg_v), "agg_v"] = 0
    summarized_hic = summarized_hic %>% filter(binA - binB != 0, agg_m >= quantile(agg_m,
                     mean_thres), agg_v >= quantile(agg_v, var_thres)) %>% select(chrom, binA, binB)
    cell_names = unique(hic_df$cell)
    input_mat = matrix(0, nrow = length(cell_names), ncol = nrow(summarized_hic))
    for (i in 1:length(cell_names)) {
      output_cell = summarized_hic
      output_cell$BandNorm = 0
      temp = hic_df %>% filter(diag > 0, cell == cell_names[i], chrom %in%
                                         summarized_hic$chrom, binA %in% summarized_hic$binA, binB %in%
                                         summarized_hic$binB)
      setDT(output_cell)
      setDT(temp)
      output_cell = output_cell[temp, `:=`(BandNorm, i.BandNorm), on = .(chrom,
                                                                         binA, binB)]
      input_mat[i, ] = output_cell$BandNorm
    }
  } else {
    paths = list.files(path, recursive = TRUE, full.names = TRUE)
    names = list.files(path, recursive = TRUE)
    cell_names = names
    load_cell = function(i) {
      return(fread(paths[i]) %>% rename(chrom = V1, binA = V2, binB = V4, BandNorm = V5) %>%
               mutate(diag = abs(binB - binA), BandNorm_s = BandNorm^2) %>% select(-V3))
    }
    summarized_hic = c()
    for (i in 1:length(path)) {
      summarized_hic = bind_rows(summarized_hic, load_cell(i)) %>% group_by(chrom,
                                                                            binA, binB, diag) %>% summarise_all(sum)
    }
    summarized_hic = summarized_hic %>% mutate(agg_v = (BandNorm_s - BandNorm^2/length(paths))/(length(paths) -
                                                                                                  1)) %>% rename(agg_m = BandNorm)
    summarized_hic = summarized_hic %>% filter(diag > 0, agg_m >= quantile(agg_m,
                                                                           mean_thres), agg_v >= quantile(agg_v, var_thres)) %>% select(chrom, binA,
                                                                                                                                        binB)
    input_mat = matrix(0, nrow = length(names), ncol = nrow(summarized_hic))
    for (i in 1:length(paths)) {
      output_cell = summarized_hic
      output_cell$BandNorm = 0
      temp = load_cell(i) %>% filter(diag > 0, chrom %in% summarized_hic$chrom,
                                     binA %in% summarized_hic$binA, binB %in% summarized_hic$binB)
      setDT(output_cell)
      setDT(temp)
      output_cell = output_cell[temp, `:=`(BandNorm, i.BandNorm), on = .(chrom,
                                                                         binA, binB)]
      input_mat[i, ] = output_cell$BandNorm
    }
  }
  pca_mat = prcomp(input_mat)$x[, 1:dim_pca]
  # Whether to use Harmony to clean the PCA embedding.
  if (do_harmony) {
    pca_mat <- HarmonyMatrix(pca_mat, batch, "dataset", do_pca = FALSE)
    pca_mat = (pca_mat)[, 1:dim_pca]
  }
  rownames(pca_mat) = cell_names
  return(pca_mat)
}

#' plot_embedding
#'
#' This function allows you to calculate and plot the UMAP or tSNE embedding.
#' @param embedding The PCA embedding obtained from create_embedding.
#' @param type A string of "tSNE" or "UMAP".
#' @param cell_type Include the cell type information in the plot. Default is NULL.
#' @export
#' @import umap
#' @import Rtsne
#' @import ggplot2
#' @examples
#' data("hic_df")
#' data("batch")
#' data("cell_type")
#' bandnorm_result = bandnorm(hic_df = hic_df, save = FALSE)
#' embedding = create_embedding(hic_df = bandnorm_result, do_harmony = TRUE, batch = batch)
#' plot_embedding(embedding, "UMAP", cell_type = cell_type)
plot_embedding = function(embedding, type, cell_type = NULL) {
  if (is.null(cell_type)){
    if (type == "tSNE") {
      embedding = Rtsne(embedding)$Y
      embedding = data.frame(embedding)
      colnames(embedding) = c("X1", "X2")
      ggplot(embedding, aes(x = X1, y = X2)) + geom_point() + ylab("tSNE 2") +
        xlab("tSNE 1") + theme_bw(base_size = 15)
    } else {
      embedding = umap(embedding)$layout
      embedding = data.frame(embedding)
      colnames(embedding) = c("X1", "X2")
      ggplot(embedding, aes(x = X1, y = X2)) + geom_point() + ylab("UMAP 2") +
        xlab("UMAP 1") + theme_bw(base_size = 15)
    }
  } else {
    if (type == "tSNE") {
      embedding = Rtsne(embedding)$Y
      embedding = data.frame(embedding)
      colnames(embedding) = c("X1", "X2")
      ggplot(embedding, aes(x = X1, y = X2)) + geom_point(aes(color = cell_type)) +
        ylab("tSNE 2") + xlab("tSNE 1") + theme_bw(base_size = 15) + labs(color = "Cell Type")
    } else {
      embedding = umap(embedding)$layout
      embedding = data.frame(embedding)
      colnames(embedding) = c("X1", "X2")
      ggplot(embedding, aes(x = X1, y = X2)) + geom_point(aes(color = cell_type)) +
        ylab("UMAP 2") + xlab("UMAP 1") + theme_bw(base_size = 15) + labs(color = "Cell Type")
    }
  }
}

