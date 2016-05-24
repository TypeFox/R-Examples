reg.dim.est.treshold = function(n,norm,Kconst=NULL){
  if (is.null(Kconst))
    Kconst = 1
  (norm*log(n))/(Kconst*n^{0.5})
}

#' Consider a linear regression problem for a multivariate stationary time series X_t:
#' \deqn{Y_t = P X_t + \varepsilon_t.}
#' Estimator based on formula
#' \deqn{EY_0X_0 (EX_0^2)^{-1} = P }
#' is fragile on the eigendirections of \eqn{EX_0^2} with small eigenvalues.
#' It is therefore desired to truncate the inversion at a level where eigenvalues are
#' estimated consistently.
#' Procedure \code{dim.est} suggest such level by taking only the eigenvalues which are greater
#' and equal than \eqn{1/\sqrt{n}K_{const}}.
#' It is designed for \eqn{\sqrt{n}} consistent matrix estimator and can serve as one of heuristics
#' for matrix inverion problems.
#' 
#' @title Estimate the optimal dimension in linear regression problem
#' @param eigenvalues vector of eigenvalues
#' @param n used for estimation
#' @param Kconst parameter for fitting the convergence rate to 1/(Kconst*n^{1/2})
#' @return number of 'safe' eigendirections
#' @seealso \code{\link{reg.est}}, \code{\link{pseudoinverse}}
#' @export
#' @references Siegfried Hormann and Lukasz Kidzinski
#' A note on estimation in Hilbertian linear models
#' Research report, 2012
#' @examples
#' n = 100
#' X = rar(n)
#' Y = rar(n)
#' CV = lagged.cov(X,Y)
#' E = eigen(CV)
#' K = reg.dim.est(E$values, n)
reg.dim.est = function(eigenvalues,n,Kconst=1){
  if (!is.numeric(n))
    stop("n must be an integer")
  if (!is.numeric(eigenvalues))
    stop("eigenvalue must be numeric")
  a = abs(eigenvalues)
  
	v = (a > reg.dim.est.treshold(n,a[1],Kconst))
	max(sum(v),1)
}