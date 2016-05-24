#' @title L2 matrix norm
#' 
#' @description Compute the L2 matrix norm of M
#' @usage l2norm(M)
#' @param M the matrix (real or complex valued)
#' 
#' @export
l2norm <- function(M) {
  s <- sqrt(spectralRadius(t(M) %*% M))
  return(s)
}

#' @title L1 matrix norm
#' 
#' @description Compute the L1 matrix norm of M
#' @usage l1norm(M)
#' @param M the matrix (real or complex valued)
#' 
#' @export
l1norm <- function(M) {
  c <- max(colSums(Mod(M)))
  return(c)
}

#' @title L-infinity matrix norm
#' 
#' @description Compute the L-infinity matrix norm of M
#' @usage lInftyNorm(M)
#' @param M the matrix (real or complex valued)
#' 
#' @export
lInftyNorm <- function(M) {
  c <- max(rowSums(Mod(M)))
  return(c)
}

#' @title Max-norm of a matrix
#' 
#' @description Compute the max-norm of M
#' @usage maxNorm(M)
#' @param M the matrix (real or complex valued)
#' 
#' @export
maxNorm <- function(M) {
  return(max(abs(M)))
}

#' @title Froebenius norm of a matrix
#' 
#' @description Compute the Froebenius norm of M
#' @usage frobNorm(M)
#' @param M the matrix (real or complex valued)
#' 
#' @export
frobNorm <- function(M) {
  A <- (t(M) %*% M)
  A <-  A * diag(nrow(A))
  return(sqrt(sum(A)))
}

#' @title Spectral radius 
#' 
#' @description Compute the spectral radius of M
#' @usage spectralRadius(M)
#' @param M the matrix (real or complex valued)
#' 
#' @export
spectralRadius <- function(M) {
  e <- eigen(M)
  maxEig <- max(Mod(e$values))
  return(maxEig)
}

#' @title Spectral norm 
#' 
#' @description Compute the spectral norm of M
#' @usage spectralNorm(M)
#' @param M the matrix (real or complex valued)
#' 
#' @export
spectralNorm <- function(M) {
  return(sqrt(spectralRadius(t(M) %*% M)))
}
