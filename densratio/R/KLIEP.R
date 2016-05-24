#' Estimate Density Ratio p_nu(x)/p_de(y) by KLIEP (Kullback-Leibler Importance Estimation Procedure)
#'
#' @param x numeric vector or matrix as data from a numerator distribution p_nu(x).
#' @param y numeric vector or matrix as data from a denominator distribution p_de(y).
#' @param sigma positive numeric vector as a search range of Gaussian kernel bandwidth.
#' @param kernel_num positive integer as number of kernels.
#' @param fold positive integer as a numer of the folds of cross validation.
#' @param verbose logical(default TRUE).
#'
#' @return KLIEP object that contains the function to estimate density ratio.
#'
#' @export
KLIEP <- function(x, y, sigma = "auto", kernel_num = 100, fold = 5, verbose = TRUE) {
  if(verbose) message("############ Start KLIEP ############")
  if(is.vector(x)) x <- matrix(x)
  if(is.vector(y)) y <- matrix(y)
  if(ncol(x) != ncol(y)) stop("x and y must be same dimensions.")

  nx <- nrow(x)

  kernel_num <- min(kernel_num, nx)
  centers <- x[sample(nx, size = kernel_num), , drop = FALSE]

  if(identical(sigma, "auto")) {
    if(verbose) message("Searching optimal sigma and lambda...")
    sigma <- KLIEP_search_sigma(x, y, centers, fold, verbose)
    if(verbose) message(sprintf("Found optimal sigma = %.5f.", sigma))
  } else if(length(sigma) != 1) {
    if(verbose) message("Searching optimal sigma and lambda...")
    sigma <- KLIEP_search_sigma_list(x, y, centers, sigma, fold, verbose)
    if(verbose) message(sprintf("Found optimal sigma = %.5f.", sigma))
  }

  if(verbose) message("Optimizing alpha...")
  phi_x <- compute_kernel_Gaussian(x, centers, sigma)
  phi_y <- compute_kernel_Gaussian(y, centers, sigma)
  alpha <- KLIEP_optimize_alpha(phi_x, phi_y)
  if(verbose) message("End.")

  result <- list(alpha = alpha,
                 kernel_info = list(
                   kernel = "Gaussian RBF",
                   kernel_num = kernel_num,
                   sigma = sigma,
                   centers = centers
                 ),
                 fold = fold,
                 compute_density_ratio = function(x) {
                   if(is.vector(x)) x <- matrix(x)
                   phi_x <- compute_kernel_Gaussian(x, centers, sigma)
                   density_ratio <- as.vector(phi_x %*% alpha)
                   density_ratio
                 }
  )
  class(result) <- c("KLIEP", class(result))
  if(verbose) message("############ Finished KLIEP #########")
  result
}
