#' Create multivariate ellipse.
#' Randomly sample points from a probability contour of a multivariate normal.
#' 
#' There are two ways to use this function.  You can either supply
#' a data set for which a multivariate normal ellipse will be drawn
#' or you can supply the mean vector, covariance matrix and number
#' of dimensions yourself. 
#'
#' @param data data frame or matrix
#' @param npoints number of points to sample
#' @param cl proportion of density contained within ellipse
#' @param mean mean vector
#' @param cov variance-covariance matrix
#' @param df degrees of freedom used for calculating F statistic
#' @export
#' @keywords internal
ellipse <- function(data, npoints=1000, cl=0.95, mean=colMeans(data), cov=var(data), df=nrow(data)) 
{
  norm.vec <- function(x) x / sqrt(sum(x^2))

  p <- length(mean)
  ev <- eigen(cov)

  sphere <- matrix(rnorm(npoints*p), ncol=p)
  cntr <- t(apply(sphere, 1, norm.vec))

  cntr <- cntr %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
  cntr <- cntr * sqrt(p * (df-1) * qf(cl, p, df-p) / (df * (df-p)))
  if (!missing(data)) colnames(cntr) <- colnames(data)

  cntr + rep(mean, each=npoints) 
}
