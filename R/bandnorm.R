#' Download existing single-cell raw data
#'
#' This function allows you to download the currently existing single cell Hi-C data. The data is 1mb resolution, and there are basic information for each of the cell, such as cell type, batch, etc.
#' @param cell_line Must be one of "Kim2020", "Lee2019", "Li2019", "Ramani2017".
#' @param cell_type If you need to download a specific cell-type from one cell line, indicate the name of the cell-type in here.
#' @param cell_path Indicate the output path for raw data.
#' @param summary_path Indicate the output path for summary data containing information of batch, cell type, depth and sparsity. Default is NULL.
#' @export
#' @examples
#' # download_schic("Li2019", cell_path = getwd())
download_schic = function(cell_line, cell_type = NULL, cell_path, summary_path = NULL) {
  if (!cell_line %in% c("Kim2020", "Lee2019", "Li2019", "Ramani2017")){
    stop("We currently don't support other cell lines. Please use one of Kim2020, Lee2019, Li2019, Ramani2017.")
  }
  if (!dir.exists(cell_path)){
    warning("path for saving the data doesn't exist, will create a path.")
    dir.create(cell_path, recursive = TRUE)
  }
  if (!is.null(summary_path)) {
    if (!dir.exists(summary_path)){
      warning("path for saving the summary data doesn't exist, will create a path.")
      dir.create(summary_path, recursive = TRUE)
    }
    input_summary = paste("http://pages.stat.wisc.edu/~sshen82/bandnorm/Summary/",
                          cell_line, "_Summary.txt", sep = "")
    invisible(download.file(input_summary, destfile = paste(summary_path, "/", cell_line,
                                                  "_Summary.txt", sep = ""), quiet = TRUE))
  }
  if (is.null(cell_type)) {
    input = paste("http://pages.stat.wisc.edu/~sshen82/bandnorm/Summary/", cell_line,
                  "_list.txt", sep = "")
    invisible(download.file(input, destfile = paste(cell_path, "/", cell_line,
                                          "_list.txt", sep = ""), quiet = TRUE))
    list_files = fread(paste(cell_path, "/", cell_line, "_list.txt", sep = ""), header = FALSE)
  } else {
    input = paste("http://pages.stat.wisc.edu/~sshen82/bandnorm/Summary/", cell_line,
                  "_", cell_type, "_list.txt", sep = "")
    invisible(download.file(input, destfile = paste(cell_path, "/", cell_line,
                                          "_", cell_type, "_list.txt", sep = ""), quiet = TRUE))
    list_files = fread(paste(cell_path, "/", cell_line, "_", cell_type, "_list.txt", sep = ""), header = FALSE)
  }
  pb = progress_bar$new(total = length(list_files$V1))
  for (i in list_files$V1){
    pb$tick()
    invisible(download.file(i, destfile = paste(cell_path, "/", basename(i), sep = ""), quiet = TRUE))
  }
  if (is.null(cell_type)) {
    invisible(file.remove(paste(cell_path, "/", cell_line, "_list.txt", sep = "")))
  } else {
    invisible(file.remove(paste(cell_path, "/", cell_line, "_", cell_type, "_list.txt", sep = "")))
  }
}

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
  if (is.null(path) && is.null(hic_df)){
    stop("Please specify one of path or hic_df.")
  }
  if (!is.null(path) && !is.null(hic_df)){
    warning("Specified both path and hic_df. Will use path file as default for saving memory.")
  }
  if (save) {
    if (is.null(save_path)){
      stop("path for saving the data is NULL!")
    }
    if (!dir.exists(save_path)){
      warning("path for saving the data doesn't exist, will create one according to the directory.")
      dir.create(save_path, recursive = TRUE)
    }
  }
  # Get path and name for all the cells in this path.
  if (!is.null(path)){
    if (!dir.exists(path)){
      stop("path for data doesn't exist!")
    }
    paths = list.files(path, recursive = TRUE, full.names = TRUE)
    names = basename(list.files(path, recursive = TRUE))
    # The input format of the cell should be [chr1, bin1, chr2, bin2, count].
    load_cell = function(i) {
      return(fread(paths[i], select = c(1, 2, 4, 5)) %>% dplyr::rename(chrom = V1, binA = V2, binB = V4, count = V5) %>%
               mutate(diag = abs(binB - binA), cell = names[i]))
    }
    hic_df = rbindlist(lapply(1:length(paths), load_cell))
  }
  # Calculate band depth and the mean of band depth for bandnorm.
  band_info <- hic_df %>% group_by(chrom, diag, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(chrom, diag) %>% summarise(depth = mean(band_depth))

  hic_df <- hic_df %>% left_join(alpha_j, by = c("chrom", "diag")) %>% left_join(band_info,
                                                                                 by = c("chrom", "diag", "cell")) %>% mutate(BandNorm = count/band_depth * depth) %>% 
  select(-c(band_depth, depth, count))
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

#' Juicer .hic version of BandNorm
#'
#' This function allows you to calculate the BandNorm normalization using the inputs that are Juicer .hic format, and the output won't be different from the original version.
#' @param path The path for all the cells in a directory. There can be sub-directories.
#' @param save Whether to save each normalized cells. Default is TRUE. Note that if don't have large memory on your computer, and you need to use create_embedding function, it is highly recommended to save the cells because it helps lower the cost of memory in this function.
#' @param save_path Indicate the output path for normalized cells. Only need it when "save" parameter is TRUE. Default is NULL.
#' @param resolution Specify the resolution from the hic file."
#' @param chroms Specify the chromosomes to use, the format is a string like "1". If the input is "all", it will include all the intra-chromosomal bin pairs. This does not allow inputting chrM and chrY currently.
#' @export
#' @import data.table
#' @import dplyr
#' @examples
#' # We will have a thorough example for this part!
bandnorm_juicer = function(path = NULL, resolution, chroms, save = TRUE, save_path = NULL) {
  if (is.null(path)){
    stop("Please specify the path.")
  }
  if (save) {
    if (is.null(save_path)){
      stop("path for saving the data is NULL!")
    }
    if (!dir.exists(save_path)){
      warning("path for saving the data doesn't exist, will create one according to the directory.")
      dir.create(save_path, recursive = TRUE)
    }
  }
  # Get path and name for all the cells in this path.
  if (!dir.exists(path)){
    stop("path for data doesn't exist!")
  }
  paths = list.files(path, recursive = TRUE, full.names = TRUE)
  if (chroms == "all"){
    chroms = strawr::readHicChroms(paths[1])$name
    chroms = chroms[!chroms %in% c("ALL", "M", "Y")]
  }
  chroms = as.character(chroms)
  names = basename(list.files(path, recursive = TRUE))
  names = gsub(".hic", ".txt", names)
  load_cell = function(i) {
    cell = list()
    cell$contact = list()
    for (j in 1:length(chroms)){
      cell$contact[[chroms[j]]] =  as.matrix(strawr::straw("NONE", paths[i], 
                                                           chroms[j], chroms[j], unit = "BP", resolution))
      colnames(cell$contact[[chroms[j]]]) = c("x", "y", "counts")
    }
    clean_cell = function(k) {
      temp = as.data.frame(cell$contact[[chroms[k]]])
      n = nrow(cell$contact[[chroms[k]]])
      temp$chromA = rep(paste("chr", chroms[k], sep = ""), n)
      return(temp %>% select(chromA, x, y, counts) %>% dplyr::rename(chrom = chromA, binA = x, binB = y, count = counts))
    }
    return(rbindlist(lapply(1:length(chroms), clean_cell)) %>%
      mutate(diag = abs(binB - binA), cell = names[i]))
  }
  hic_df = rbindlist(lapply(1:length(paths), load_cell))
  # Calculate band depth and the mean of band depth for bandnorm.
  band_info <- hic_df %>% group_by(chrom, diag, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(chrom, diag) %>% summarise(depth = mean(band_depth))

  hic_df <- hic_df %>% left_join(alpha_j, by = c("chrom", "diag")) %>% left_join(band_info,
                                                                                 by = c("chrom", "diag", "cell")) %>% mutate(BandNorm = count/band_depth * depth) %>%
  select(-c(band_depth, depth, count))
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

#' Create PCA embedding
#'
#' This function allows you to obtain the PCA embedding for tSNE or UMAP plotting.
#' @param path The path for all the normalized cells in a directory. There can be sub-directories. Using path means loading the bandnorm normalized data iteratively from the directory, so it is relatively slower than using hic_df. However, it won't eat up your memory too much, and the speed is acceptable.
#' @param hic_df After using bandnorm, if you keep the data frame, it is possible to use this instead of "path" as input of create_embedding. If using hic_df, the speed will be faster, but is costs more memory.
#' @param chrs Chromosomes used in the embedding. Default is chr1, ..., chr19, and chrX.
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
create_embedding = function(path = NULL, hic_df = NULL, chrs = paste0("chr", c(1:19, "X")),
                            dim_pca = 50, do_harmony = FALSE, batch = NULL) {
  # Function to create embedding for cells, after combining all the bin-pairs from all chromosomes,
  # we do PCA first, and this function returns the PCA embedding.
  # There are two options:
  # If using hic_df, the speed will be faster, but is costs more memory, so you can use it if the memory
  # of the computer is sufficient;
  # For path, it means loading the bandnorm normalized data iteratively, so it is slower. However,
  # it won't eat up your memory too much, and the speed is also acceptable.
  if (is.null(path) && is.null(hic_df)){
    stop("Please specify one of path or hic_df.")
  }
  if (!is.null(path) && !is.null(hic_df)){
    warning("Specified both path and hic_df. Will use hic_df file as default.")
  }
  if (is.null(path)) {
    n = length(unique(hic_df$cell))
    if (dim_pca > n){
      stop("The dimension for the PCA embedding should be smaller than the size of the dataset.")
    }
    if (!is.null(batch)){
      if (nrow(batch) != n){
        stop("The number of cells in batch file is not the same with the cell in path or hic_df.")
      }
      if (ncol(batch) != 2){
        stop("The batch file should only contain two columns, the first is cell names and the second is batch information.")
      }
    }
    setDT(hic_df)
    hic_df = hic_df[hic_df$diag != 0, ]
    names = unique(hic_df$cell)
    n = length(names)
    hic_df$cell = as.numeric(factor(hic_df$cell, level = names))
    hic_df = hic_df %>% 
      mutate(featureIndex = paste(binA, binB, sep = "_")) %>%
      select(-c(binA, binB, diag))
    output = list()
    for (c in 1:length(chrs)) {
      hic_dfCHR = hic_df[chrom == chrs[c], ]
      hic_dfCHR$featureIndex = as.numeric(factor(hic_dfCHR$featureIndex))
      hic_dfCHR = sparseMatrix(i = hic_dfCHR$cell, j = hic_dfCHR$featureIndex, x = hic_dfCHR$BandNorm, 
                               dims = c(max(hic_df$cell), max(hic_dfCHR$featureIndex)),
                               index1 = TRUE)
      output[[chrs[c]]] = hic_dfCHR
    }
  } else {
    paths = list.files(path, recursive = TRUE, full.names = TRUE)
    names = basename(list.files(path, recursive = TRUE))
    n = length(names)
    if (dim_pca > n){
      stop("The dimension for the PCA embedding should be smaller than the size of the dataset.")
    }
    if (!is.null(batch)){
      if (nrow(batch) != n){
        stop("The number of cells in batch file is not the same with the cell in path or hic_df.")
      }
      if (ncol(batch) != 2){
        stop("The batch file should only contain two columns, the first is cell names and the second is batch information.")
      }
    }
    load_cell = function(i) {
      return(fread(paths[i], select = c(1, 2, 4, 5))[V2 - V4 != 0 & V1 == chrs[c]] %>% 
               mutate(cellIndex = i, featureIndex = paste(V2, V4, sep = "_")) %>%
               select(-c(V1, V2, V4)))
    }
    output = list()
    for (c in 1:length(chrs)) {
      hic_df = rbindlist(lapply(1:length(paths), load_cell))
      hic_df$featureIndex = as.numeric(as.factor(hic_df$featureIndex))
      hic_df = sparseMatrix(i = hic_df$cellIndex, j = hic_df$featureIndex, x = hic_df$V5, 
                            dims = c(max(hic_df$cellIndex), max(hic_df$featureIndex)),
                            index1 = TRUE)
      output[[chrs[c]]] = hic_df
    }
  }
  if (dim(output[["chr1"]])[2] < 80000){
    SVDInput = c()
    for (i in 1:length(output)){
      SVDInput = cbind(SVDInput, output[[i]])
    }
    pca_mat = sparsesvd(SVDInput, dim_pca)
    pca_mat = pca_mat$u %*% diag(pca_mat$d)
  } else {
    pca_mat = c()
    for (i in 1:length(output)){
      SVDInput = sparsesvd(output[[i]], dim_pca)
      SVDInput = SVDInput$u %*% diag(SVDInput$d)
      pca_mat = cbind(pca_mat, SVDInput)
    }
    pca_mat = fast.prcomp(pca_mat)$x[, 1:dim_pca]
  }
  # Whether to use Harmony to clean the PCA embedding.
  if (do_harmony) {
    colnames(batch) = c("cell_name", "batch")
    batch = data.frame(batch)
    batch = batch[match(names, batch$cell_name), ]$batch
    pca_mat <- harmony::HarmonyMatrix(pca_mat, batch, "dataset", do_pca = FALSE)
    pca_mat = (pca_mat)[, 1:dim_pca]
  }
  rownames(pca_mat) = names
  return(pca_mat)
}

#' Plot tSNE or UMAP embedding
#'
#' This function allows you to calculate and plot the UMAP or tSNE embedding.
#' @param embedding The PCA embedding obtained from create_embedding.
#' @param type A string of "tSNE" or "UMAP". Default is "UMAP".
#' @param cell_info Include the cell information in the plot. It should be a matrix or data.frame, and the first column is cell name in the corresponding to the embedding, and the second column is the information. Default is NULL.
#' @param label A string indicating the name of the cell information you want to include. Default is NULL.
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
#' plot_embedding(embedding, "UMAP", cell_info = cell_type, label = "Cell Type")
plot_embedding = function(embedding, type = "UMAP", cell_info = NULL, label = NULL) {
  if (!type %in% c("UMAP", "tSNE")){
    stop("We currently only support for UMAP or tSNE embedding.")
  }
  if (is.null(cell_info)){
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
    if (nrow(cell_info) != nrow(embedding)){
      stop("The number of cells in batch file is not the same with the number of cells the embedding.")
    }
    if (ncol(cell_info) != 2){
      stop("The batch file should only contain two columns, the first is cell names and the second is batch information.")
    }
    if (is.null(label)){
      warning("The label is empty. Will substitute with the default 'cell information'")
      label = "cell information"
    }
    cell_names = rownames(embedding)
    colnames(cell_info) = c("cell_name", "info")
    cell_info = data.frame(cell_info)
    cell_info = cell_info[match(cell_names, cell_info$cell_name), ]$info
    if (type == "tSNE") {
      embedding = Rtsne(embedding)$Y
      embedding = data.frame(embedding)
      colnames(embedding) = c("X1", "X2")
      ggplot(embedding, aes(x = X1, y = X2)) + geom_point(aes(color = cell_info)) +
        ylab("tSNE 2") + xlab("tSNE 1") + theme_bw(base_size = 15) + labs(color = label)
    } else {
      embedding = umap(embedding)$layout
      embedding = data.frame(embedding)
      colnames(embedding) = c("X1", "X2")
      ggplot(embedding, aes(x = X1, y = X2)) + geom_point(aes(color = cell_info)) +
        ylab("UMAP 2") + xlab("UMAP 1") + theme_bw(base_size = 15) + labs(color = label)
    }
  }
}


#' Cooler hic version of BandNorm (Thanks skelviper for supporting the code!)
#'
#' This function allows you to calculate the BandNorm normalization using the inputs that are cooler .mcool format. (We might also add .scool format.)
#' @examples
#' # Will add example here.


bandnorm_cooler <- function(coolerPath, resolution, cellname = NULL){
  if(length(cellname)==0){
    cellname = str_split(coolerPath,"/")[[1]] %>% tail(n=1) %>% str_remove(".mcool")
  }
  if(str_split(coolerPath,fixed("."))[[1]] %>% tail(n=1) == "mcool"){
    dump <- rhdf5::h5dump(coolerPath)
    dump <- dump$resolutions[[which(names(dump$resolutions) ==format(resolution, scientific = FALSE))]]
  }
  
  ids <- data.frame(chr = dump$bins$chrom, start = dump$bins$start, end = dump$bins$end, id = seq(1, length(dump$bins$chrom), by = 1)-1)
  mat <- data.frame(bin1 = dump$pixels$bin1_id, bin2 = dump$pixels$bin2_id, IF = dump$pixels$count)
  colnames(mat) <- c('i', 'j', 'IF')
  colnames(ids) <- c('chr1', 'start1', 'end1', 'id')
  
  new_mat <- left_join(mat, ids, by = c('i' = 'id'))
  colnames(ids) <- c('chr2', 'start2', 'end2', 'id')
  new_mat <- left_join(new_mat, ids, by = c('j' = 'id'))
  new_mat <- new_mat %>% filter(chr1==chr2)
  new_mat <- new_mat %>% mutate(diag = (j-i)*resolution) %>% select(4,5,8,3,diag) %>% mutate(cellname=cellname)
  names(new_mat) <- c("chrom","binA","binB","count","diag","cell")
  return(new_mat)
}

