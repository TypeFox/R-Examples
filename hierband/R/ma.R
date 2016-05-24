
#' Covariance of an equal-weighted moving-average process
#'
#' Here, \code{Sig[j,k] = 0} if \code{|j-k|>K}
#' and \code{Sig[j,k] = 1 - |j-k| / K} otherwise.
#'
#' @param p dimension of covariance matrix
#' @param K moving-average bandwidth
#' 
#' @return Returns the covariance matrix, \code{Sig}, and the symmetric
#' square root, \code{A}, of this matrix.
#'
#' @export
ma <- function(p, K) {
  sig <- seq(K,1) / K
  banded(p, sig)
}

#' Generates a banded covariance matrix and matrix squareroot
#' sig: value of kth band (starting with diagonal)
#' size of band is length(sig)
#' 
#' @param p dimension of covariance matrix
#' @param sig vector of values of Toeplitz matrix
banded <- function(p, sig) {
  ii <- toeplitz(1:p)
  K <- length(sig)
  Sig <- ii <= K
  for (k in seq(K)) Sig[ii==k] <- sig[k]

  eig <- eigen(Sig)
  A <- eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)
  list(Sig=Sig, A=A)
}
