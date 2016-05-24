geigen <- function(Amat, Bmat, Cmat)
{
  #  solve the generalized eigenanalysis problem
  #
  #    max {tr L'AM / sqrt[tr L'BL tr M'CM] w.r.t. L and M
  #
  #  Arguments:
  #  AMAT ... p by q matrix
  #  BMAT ... order p symmetric positive definite matrix
  #  CMAT ... order q symmetric positive definite matrix
  #  Returns:
  #  VALUES ... vector of length s = min(p,q) of eigenvalues
  #  LMAT   ... p by s matrix L
  #  MMAT   ... q by s matrix M

  #  last modified 2007 to use svd2;  previous mod 18 May 2001

  Bdim <- dim(Bmat)
  Cdim <- dim(Cmat)
  if (Bdim[1] != Bdim[2]) stop('BMAT is not square')
  if (Cdim[1] != Cdim[2]) stop('CMAT is not square')
  p <- Bdim[1]
  q <- Cdim[1]
  s <- min(c(p,q))
  if (max(abs(Bmat - t(Bmat)))/max(abs(Bmat)) > 1e-10) stop(
    'BMAT not symmetric.')
  if (max(abs(Cmat - t(Cmat)))/max(abs(Cmat)) > 1e-10) stop(
    'CMAT not symmetric.')
  Bmat  <- (Bmat + t(Bmat))/2
  Cmat  <- (Cmat + t(Cmat))/2
  Bfac  <- chol(Bmat)
  Cfac  <- chol(Cmat)
  Bfacinv <- solve(Bfac)
  Cfacinv <- solve(Cfac)
  Dmat <- t(Bfacinv) %*% Amat %*% Cfacinv
  if (p >= q) {
    result <- svd2(Dmat)
    values <- result$d
    Lmat <- Bfacinv %*% result$u
    Mmat <- Cfacinv %*% result$v
  } else {
    result <- svd2(t(Dmat))
    values <- result$d
    Lmat <- Bfacinv %*% result$v
    Mmat <- Cfacinv %*% result$u
  }
  geigenlist <- list (values, Lmat, Mmat)
  names(geigenlist) <- c('values', 'Lmat', 'Mmat')
  return(geigenlist)
}
