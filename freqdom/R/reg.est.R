#' Estimate regression operator \eqn{P} (matrix \eqn{d \times d}) in model 
#' \deqn{Y_t = P X_t + \varepsilon_t,} where \eqn{X_t} is a \eqn{d}-dimensional stationary process
#' and \eqn{\varepsilon_t} forms a white noise.
#'
#' @title Estimate the optimal dimension in linear regression problem
#' @param X regressors process
#' @param Y response process
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
#' Psi = reg.est(X,Y)
#' @seealso \code{\link{reg.dim.est}}, \code{\link{speclagreg}}
reg.est = function (X,Y,K=NULL,Kconst=1){
  if (is.vector(X))
    X = matrix(X)
  if (is.vector(Y))
    Y = matrix(Y)
  
	nbasis = dim(Y)[1]
	n = dim(Y)[2]

	Cv = (t(X) %*% X) / n
	Dv = (t(Y) %*% X) / n
	E = eigen(Cv)

	if (is.null(K))
		th = reg.dim.est.treshold(n,E$values[1],Kconst)
  else
    th = E$values[K] - 0.00001

	Dv %*% pseudoinverse(Cv,sum(abs(E$values) >= th))
}
