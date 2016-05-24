#' @title Banding Opreator on Covariance Matrix
#'
#' @description
#' Apply banding operator on a covariance matrix with a banding parameter.
#'
#' @param sigma a p*p covariance matrix
#' @param k banding parameter
#' @return a regularized covariance matrix after banding operation
#' @references "High-Dimensional Covariance Estimation" by Mohsen Pourahmadi
#' @examples
#' data(m.excess.c10sp9003)
#' cov.SAM <- cov(m.excess.c10sp9003)
#' banding(cov.SAM, 7)
#' @export

banding <- function(sigma, k = 0) {
  p <- ncol(sigma)
  for (i in 1:p) {
    for (j in 1:i) {
      if (abs(i-j) > k) sigma[i,j] <- 0
      sigma[j,i] <- sigma[i,j]
    }
  }
  return(sigma)
}