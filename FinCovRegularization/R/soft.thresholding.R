#' @title Soft-Thresholding Opreator on Covariance Matrix
#'
#' @description
#' Apply soft-thresholding operator on a covariance matrix with 
#' a soft-thresholding parameter.
#'
#' @param sigma a covariance matrix
#' @param threshold soft-thresholding parameter
#' @return a regularized covariance matrix after soft-thresholding operation
#' @references "High-Dimensional Covariance Estimation" by Mohsen Pourahmadi
#' @examples
#' data(m.excess.c10sp9003)
#' cov.SAM <- cov(m.excess.c10sp9003)
#' soft.thresholding(cov.SAM, threshold = 0.001)
#' @export

soft.thresholding <- function(sigma, threshold = 0.5) {
  sigma.threshold <- sigma * 0
  sigma.threshold <- sigma.threshold + diag(diag(sigma))
  sigma <- sigma - diag(diag(sigma.threshold))
  sigma.threshold <- sigma.threshold + sign(sigma) * pmax(abs(sigma) - threshold, 0)
  return(sigma.threshold)
}