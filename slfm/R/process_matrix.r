#' Pre-process data for SLFM
#'
#' This function pre-process the data to be used
#' for fitting a sparse latent factor model.
#'
#' @param path path containing the set of matrices to be processed
#' @param output_path path to save the processed matrices
#' @param sample_size number of matrices to be used on the principal component analysis
#' @importFrom stats cov
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export

process_matrix <- function(path, output_path, sample_size) {

  files_list <- list.files(path, recursive = TRUE)
  if(!file.exists(output_path)) {
    dir.create(output_path)
  }
  
  sproc <- process_sample(path, files_list, sample_size)

  lapply(files_list, function(x) {
    mat <- read.table(paste(path, x, sep = ""))
    x_star <- t(apply(mat, 1, function(y) y - sproc$lcm))
    x_star <- scale(x_star)
    res <- t(apply(x_star, 1, function(y) { y - y%*%sproc$pca%*%sproc$pca }))
    file_dir <- dirname(paste(output_path,x,sep="/"))
    if(!file.exists(file_dir)) {
      dir.create(dirname(paste(output_path,x,sep="/")))
    }
    write.table(res, file=paste(output_path, x, sep="/"), row.names = F, col.names = F)
  })

  message(paste("Processed data was saved to \"", output_path, "\".", sep = ""))
}

process_sample <- function(path, files_list, sample_size) {

  samp <- sample(files_list, sample_size)

  mats <- lapply(samp, function(x) read.table(paste(path, x, sep="")))

  big_mat <- do.call(rbind, mats)
  big_mat_star <- exp(big_mat)
  
  column_means <- apply(big_mat_star, 2, mean)
  w <- t(apply(big_mat_star, 1, function(x) x/column_means))
  lcm <- log(column_means)
  lw <- log(w)
  lw <- scale(lw)

  pca <- eigen(cov(lw))$vectors[,1]
  list(pca = pca, lcm = lcm)
}