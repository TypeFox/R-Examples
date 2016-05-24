#' Estimate regresson operators in a lagged linear model using spectral methods.
#' Assume model
#' \deqn{Y_t = \sum_{k=-q}^p A_k X_{t-k} + \varepsilon_t}
#' where \eqn{X_t} is a stationary multivariate time series, \eqn{(A_k)_{-q \leq k \leq p}} is a filter and \eqn{\varepsilon_t} is white noise.
#' Function \code{speclagreg} estimates parameters \eqn{A_k} with \eqn{k \in }\code{lags}
#'
#' @title Estimate regresson operators in a lagged linear model
#' @param X first process 
#' @param Y second process, if null then autocovariance of X is computed
#' @param Kconst used for heuristic as in \code{\link{reg.dim.est}}
#' @param K dimension for inversion if no heuristic should be used
#' @param lags which \eqn{A_k} should be estimated
#' @param freq grid of frequencies for computation as in \code{\link{fourier.transform}}
#' @param p window size for estimation of spectral density of X
#' @param q window size for estimation of spectral density of Y and X
#' @param weights as in \code{\link{spectral.density}}
#' @return \code{timedom} operators
#' @seealso \code{\link{linproc}}
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' #estimate regressors in model $Y_t = \sum_{i\in Z} A_i X_{t-i}$
#' A = speclagreg(X,Y) 
#' # check an advanced examples in demo(lagged.reg)
#' @export
speclagreg = function(X,Y,Kconst=1,K=NULL,lags=0:0,freq=NULL,p=10,q=10,weights="Bartlett")
{
  if (!is.matrix(X) || !is.matrix(Y))
    stop("X and Y must be matrices")
  if (dim(X)[1] != dim(Y)[1])
    stop("Number of observations must match")

  thetas = freq
  
  if (is.null(thetas))
  {
    T = 100
    thetas = pi*((-T):T)/T
  }	

  SXX = spectral.density(X,freq=thetas,q=p,weights=weights)
  SYX = spectral.density(Y,X,freq=thetas,q=q,weights=weights)
  
  if (is.null(K)){
#    K = speclagreg.K(X,Y,freq=thetas,lags=lags,p=p,q=q,weights=weights)
    K = speclagreg.K.experimental(X,Y,freq=thetas,method=0,SXX=SXX)
  }

  R = freqdom.ratio(SYX,SXX,dim(X)[1],Kconst=Kconst,K=K) # frequencywise overloaded operator

  A = invfourier(rev(R),lags=lags)
  A$estdim = K
  A
}
