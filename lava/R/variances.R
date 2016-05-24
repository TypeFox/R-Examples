### Return position of variance elements in the parameter vector (without mean parameters)
### Optimization constraints are needed on these parameters
##' @export
variances <- function(x,mean=FALSE) {
##  if (is.null(x$parpos))
##    x$parpos <- parpos(x)
  x$parpos <- parpos(Model(x),mean=TRUE)
  res <- diag(x$parpos$P)[which(diag(index(x)$P0)==1)]
  if (!mean) {
    return(res - index(x)$npar.mean)
  }
  return(res)
}
## And the off-diagonal (covariance) parameters
##' @export
offdiags <- function(x,mean=FALSE) {
  parpos <- parpos(x,mean=mean)
  pp <- parpos$P
  pp[lower.tri(pp)][(index(x)$P0)[lower.tri(pp)]==1]
}
