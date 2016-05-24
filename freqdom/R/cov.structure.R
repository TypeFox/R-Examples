#' Given two multivariate time series \eqn{X_t} and \eqn{Y_t}
#' estimates the covariances of \eqn{Cov(X_{k} Y_0)} for \eqn{k \in [-q,q]}.
#' \code{\link{lagged.cov}} is used for the estimation at each lag.
#'
#' @title Estimate the covariance structure within a given window \eqn{k \in [-q,q]}
#' @param X first process 
#' @param Y second process, if null then autocovariance of \code{X} is computed
#' @param q size of the window (covariances from \code{-q} to \code{q} will be computed)
#' @return a time domain operator
#' @export
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' cov.structure(X,Y)
cov.structure = function(X,Y=NULL,q=10){
	if (is.null(Y))
		Y = X
  
	if (!is.matrix(X) || !is.matrix(Y))
	  stop("X and Y must be matrices")
	if (dim(X)[1] != dim(Y)[1])
	  stop("Number of observations must be equal")
	
	nbasisX = dim(X)[2]
	nbasisY = dim(Y)[2]
	n = dim(X)[1]

	Ch = array(0,c(2*q+1,nbasisX,nbasisY))
	
	for (h in (-q):q)
		Ch[h+q+1,,] = lagged.cov(X,Y,h)
	Ch
  
  timedom(Ch,-q:q)
}
