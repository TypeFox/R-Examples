#' Fit SLFM to the matrices inside a directory
#'
#' This function is used to fit a Bayesian sparse
#' latent factor model to a directory of numeric matrices.
#'
#' @param path path to read the matrices from
#' @param recursive if the function should look recursively inside folders
#' @param a prior shape parameter for Gamma distribution 
#' @param b prior scale parameter for Gamma distribution
#' @param gamma_a prior parameter for Beta distribution
#' @param gamma_b prior parameter for Beta distribution
#' @param omega_0 prior variance of the spike component
#' @param omega_1 prior variance of the slab component
#' @param sample sample size after burn-in
#' @param burnin burn-in size
#' @param lag lag for MCMC
#' @param degenerate use the degenerate version of mixture
#' @importFrom coda HPDinterval
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.table
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @export
slfm_list <- function(
  path = ".", recursive = TRUE,
  a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega_0 = 0.01, omega_1 = 10, sample = 1000, burnin = round(0.25*sample), lag = 1, degenerate = FALSE) {

  files_list <- list.files(path, recursive = recursive, full.names = T)

  message("* |Press Ctrl + C to cancel...")
  pb <- txtProgressBar(0, length(files_list), style = 3, title = "")

  MATRIX_CLASSIFICATION <- c("Present", "Marginal", "Absent")
  names(MATRIX_CLASSIFICATION) <- c("S", "I", "N")

  results_list <- list()
  for(i in 1:length(files_list)) {
    file_name <- files_list[i]
    mat <- read.table(file_name)

    res <- slfm(mat, a, b, gamma_a, gamma_b, omega_0, omega_1, sample, burnin, lag, degenerate)
    clas_table <- table(res$classification)
    final_clas <- MATRIX_CLASSIFICATION[names(which.max(clas_table))]
    freq <- format(round(clas_table["S"]/sum(clas_table), 4), nsmall = 4)
    results_list[[i]] <- c(name = basename(tools::file_path_sans_ext(file_name)), clas = final_clas, frequency = freq)

    setTxtProgressBar(pb, i)
  }
  close(pb)
  ret <- as.data.frame(do.call(rbind, results_list))
  names(ret) <- c("File", "Classification", "Significant %")
  ret
}