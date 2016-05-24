#' @title Tapering Opreator on Covariance Matrix
#'
#' @description
#' Apply tapering operator on a covariance matrix with tapering parameters.
#'
#' @param sigma a p*p covariance matrix
#' @param l tapering parameter
#' @param h the ratio between taper l_h and parameter l
#' @return a regularized covariance matrix after tapering operation
#' @references "High-Dimensional Covariance Estimation" by Mohsen Pourahmadi
#' @examples
#' data(m.excess.c10sp9003)
#' cov.SAM <- cov(m.excess.c10sp9003)
#' tapering(cov.SAM, l=7, h = 1/2)
#' @export

tapering <- function(sigma, l, h = 1/2) {
  p <- ncol(sigma)
  operator <- sigma * 0
  for (i in 1:p) {
    for (j in 1:i) {
      if (abs(i-j) <= l*h) {
        operator[i,j] <- 1
      } else if ((l*h < abs(i-j)) & (abs(i-j)<l)) {
        operator[i,j] <- 2 - (abs(i-j) / (l*h)) 
      } else {
        operator[i,j] <- 0
      }
      operator[j,i] <- operator[i,j]
    }
  }
  sigma.tapering <- sigma * operator
  return(sigma.tapering)
}