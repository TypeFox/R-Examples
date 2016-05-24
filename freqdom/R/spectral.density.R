#' Computes the spectral density of processes \eqn{X_t} and \eqn{Y_t} using
#' a Bartlett style estimator, i.e.
#' \deqn{ \hat F_\theta^{XY} = \sum_{k=-q}^q W(|k|/q) \hat C_{XY}^k e^{-i\theta k},}
#' where \eqn{\theta \in [-\pi,\pi]} and \eqn{\hat C_{XY}^k} is the estimated covariance with lag \eqn{k}. 
#' Quality of the estimation depends on choise of the window size \eqn{q}
#' and \eqn{W} (\code{weights}).
#' For details on spectral density estimation please refer to "Time Series: Theory and Methods"
#' by Peter J. Brockwell and Richard A. Davis.
#' Note that estimator is calculated on the finite grid \code{thetas} so #' in some cases
#' numerical quality can be improved by choosing a more dense set.
#'
#' @title Compute the cross spectral density of processes X and Y 
#' @param X first process 
#' @param Y second process, if \code{NULL} then spectral density of X is computed
#' @param V correlation structure between coefficients of vectors (default diagonal)
#' @param freq evaluation grid - vector of values between \code{[-pi,pi]}
#' @param q size of the window (covariances from -q to q will be computed)
#' @param weights kernel used to decay significance of covariances with higher lags ('Bartlett', 'trunc', 'Tukey', 'Parzen', 'Bohman', 'Daniell', 'ParzenCogburnDavis').
#' @return Frequency Domain Operator object
#' @references Peter J. Brockwell and Richard A. Davis
#' Time Series: Theory and Methods
#' Springer Series in Statistics, 2009
#' @export
#' @keywords spec
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' spectral.density(X,Y)
spectral.density = function(X,Y=NULL,V=NULL,freq=NULL,q=NULL,weights=NULL){
  if (is.null(Y))
    Y = X
  if (is.null(q))
    q = 10
  if (is.null(V))
    V = diag(dim(X)[2])

  if (!is.matrix(X) || !is.matrix(Y))
    stop("X and Y must be matrices")
  if (dim(X)[1] != dim(Y)[1])
    stop("Number of observations must be equal")
  if (!is.positiveint(q))
    stop("q must be a positive integer")
  
  thetas = freq
  
	nbasisX = dim(X)[2]
	nbasisY = dim(Y)[2]
	n = dim(X)[1]
	Ch = cov.structure(X,Y,q)
  for (i in 1:(q*2+1))
    Ch$operators[i,,] = Ch$operators[i,,] %*% V
	
  wfunc = weights.Bartlett
  if (is.null(weights))
	  wfunc = weights.Bartlett
  else if (weights=="Bartlett")
    wfunc = weights.Bartlett
  else if (weights=="trunc")
    wfunc = weights.trunc
  else if (weights=="Tukey")
    wfunc = weights.Tukey
  else if (weights=="Parzen")
    wfunc = weights.Parzen
  else if (weights=="Bohman")
    wfunc = weights.Bohman
  else if (weights=="Daniell")
    wfunc = weights.Daniell
  else if (weights=="ParzenCogburnDavis")
    wfunc = weights.ParzenCogburnDavis
  else
    stop(paste("No weight function called",weights))
  
  weights = wfunc(-q:q/q)
  
  for (i in 1:dim(Ch$operators)[1])
    Ch$operators[i,,] = weights[i] * Ch$operators[i,,]
  
  fourier.transform(Ch, freq=thetas)
}