tc = function(X){t(Conj(X))}

#' Determine the subspace dimension for the frequency-wise operator inversion in the
#' estimation of regression coefficients
#'
#' @title Determine the subspace dimension for inversion in \code{\link{speclagreg}}
#' @param X first process 
#' @param Y second process, if null then autocovariance of X is computed
#' @param lags which \eqn{A_k} should be estimated
#' @param freq grid of frequencies for computation as in \code{\link{fourier.transform}}
#' @param p window size for estimation of spectral density of X
#' @param q window size for estimation of spectral density of Y and X
#' @param weights as in \code{\link{spectral.density}}
#' @return Heursiticaly determined subspace size
#' @seealso \code{\link{linproc}}
#' @examples
#' n = 100
#' d = 5
#' X = rar(n,d=d,Psi=matrix(0,d,d))  			# independent d-dim variables
#' w = 0.4
#' Y = w*X + (1-w)*rar(n,d=d,Psi=matrix(0,d,d))	# independent d-dim variables
#' K = speclagreg.K(X, Y, lags=-2:2)
#' @export 
speclagreg.K = function(X,Y,lags=0:0,freq=NULL,p=10,q=10,weights="Bartlett"){  
  n = dim(X)[1]
  d = dim(X)[2]
  tr = 1:(floor(n * 0.8))
  ts = (length(tr)+1):(n)

  thetas = freq
  
  if (is.null(thetas))
  {
    T=100
    thetas = pi*((-T):T)/T
  }
  
  RES = rep(0,d+1)

  for (k in 0:d)
  {
    A = speclagreg(X[tr,],as.matrix(Y[tr,]),freq=freq,K=k,lags=lags,p=p,q=q,weights=weights)
    Yhat = A %c% X[ts,]
    RES[k+1] = sum((Yhat - as.matrix(Y[ts,]))^2)/n
  }
  
  K = which.min(RES) - 1
  K
}
