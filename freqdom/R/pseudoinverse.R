#' Invert first \code{K} eigendirections of the matrix. If \code{K} is
#' not specified the functions takes direction with eigencalues grater
#' than given treshold \code{th}.
#' If \code{th} is also not specified then all direction are inverted
#' (equivalent to \code{\link[base]{solve}})
#'
#' @title Invert first K eigendirections of the matrix.
#' @param M matrix to solve
#' @param K number of directions to invert
#' @param th fixed treshold for eigenvalues
#' @return Pseudoinverse matrix
#' @export 
pseudoinverse = function(M,K=NULL,th=NULL){
  if (!is.matrix(M) || dim(M)[1] != dim(M)[2])
    stop("M must be a square matrix")
  if (!is.null(K) && !is.positiveint(K+1))
    stop("K must be a nonnegative integer")
  
	nbasis = dim(M)[1]
	if (is.null(nbasis) || nbasis ==1) 
		res = 1 / M
	else {
		E = eigen(M)
		vals = 1 / E$values
    
		if (is.null(K))
		  K = nbasis
    if (!is.null(th))
      K = min(K, sum(abs(E$values) > th))
    
		vals[ 1:nbasis > K ] = 0

		res = E$vectors %*% diag(vals) %*% t(Conj(E$vectors))
	}
	res
}
