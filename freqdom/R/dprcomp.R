#' For a given process \code{X} eigendecompose it's spectral density
#' and use an inverse fourier transform to get coefficients of the optimal
#' filter. For details please refer to Hormann et al paper.
#'
#' @title Compute DPCA filter coefficients
#' @param X multivariate stationary time series
#' @param V correlation structure between coefficients of vectors (default diagonal)
#' @param lags requested filter coefficients
#' @param q window for spectral density estimation as in \code{\link{spectral.density}}
#' @param weights as in \code{\link{spectral.density}}
#' @param freq frequency grid to estimate on as in \code{\link{spectral.density}}
#' @return principal components series
#' @references Siegfried Hormann, Lukasz Kidzinski and Marc Hallin
#' Dynamic Functional Principal Component
#' Research report, 2012
#' @seealso \code{\link{dpca.inverse}}, \code{\link{dpca.scores}}
#' @export
dprcomp = function(X,V=NULL,lags=-10:10,q=NULL,weights=NULL,freq=NULL){
  if (!is.matrix(X))
    stop("X must be a matrix")
  if (!is.vector(lags) || any(!is.positiveint(abs(lags))))
    stop("lags must be a vector of integers")
  if (is.null(V))
    V = diag(dim(X)[2])

  SD = spectral.density(X,q=q,weights=weights,freq=freq)
  E = freqdom.eigen(SD)
  
	nbasis = dim(E$vectors)[2]

  XI = array(0,c(length(lags),nbasis,nbasis))

	for (component in 1:nbasis)
    XI[,component,] = t(exp(-(SD$freq %*% t(lags)) * 1i)) %*% E$vectors[,,component] / length(SD$freq) 
  # TODO: this is nothing but 'invfourier' written in a fancy way - should be simplified

	timedom(Re(XI[length(lags):1,,]),lags)
}
