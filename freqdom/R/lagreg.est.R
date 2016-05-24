#' Estimate regression operator \eqn{P_k} (matrices \eqn{d \times d}) in model 
#' \deqn{Y_t = \sum_{l \in L} P_k X_{t-l} + \varepsilon_t,} where \eqn{X_t} is a \eqn{d}-dimensional stationary process
#' and \eqn{\varepsilon_t} forms a white noise.
#'
#' @title Estimate the optimal dimension in linear regression problem
#' @param X regressors process
#' @param Y response process
#' @param lags lags which should be estimated
#' @param K how many directions should be inverted (as in \code{\link{pseudoinverse}})
#' @param Kconst constant for heuristic (as in \code{\link{reg.dim.est}})
#' @return Estimated regression operator
#' @export
#' @references Siegfried Hormann and Lukasz Kidzinski
#' A note on estimation in Hilbertian linear models
#' Research report, 2012
#' @examples
#' X = rar(100)
#' e = rar(100)
#' Y = X + 0.3 * e
#' Psi = lagreg.est(X,Y,lags=0:2)
#' @seealso \code{\link{reg.dim.est}}, \code{\link{speclagreg}}
lagreg.est = function (X,Y,lags=-5:5,K=NULL,Kconst=1){
  if (is.vector(X))
    X = as.matrix(X)
  if (is.vector(Y))
    Y = as.matrix(Y)
  
  nbasisX = dim(X)[2]
  n = dim(X)[1]
  nbasisY = dim(Y)[2]

  R = c()
  for (lag in lags){
    Xshifted = X
    if (lag > 0)
      Xshifted = rbind(matrix(0,lag,nbasisX),as.matrix(X[1:(n-lag),]))
    if (lag < 0){
      Xshifted = rbind(as.matrix(X[abs(lag) + 1:(n-abs(lag)),]),matrix(0,abs(lag),nbasisX))
    }
    R = cbind(R,Xshifted)
  }
  X = R

	Cv = (t(X) %*% X) / n
	Dv = (t(Y) %*% X) / n
	E = eigen(Cv)
  
  if (is.null(K))
    K = nbasisY*length(lags)
  
  RES = Dv %*% pseudoinverse(Cv,K)

  A = array(0,c(length(lags),nbasisY,nbasisX))
  for (i in 1:length(lags)){
    A[i,,] = RES[,(i-1) * nbasisX + 1:nbasisX]
  }
  
  timedom(A,lags)
}
