# We compute the solution path of the trend filtering problem:
#
# \hat{\beta}(\lambda) =
# \argmin_\beta \|y - X \beta|_2^2 + \lambda\|D \beta\|_1,
#
# where D is (p-1) x p is the discrete difference operator (of
# order 1), and X is n x p with full column rank. The solution is
# piecewise constant, with adaptively chosen break points.

fusedlasso1d <- function(y, pos, X, gamma=0, approx=FALSE, maxsteps=2000,
                         minlam=0, rtol=1e-7, btol=1e-7, eps=1e-4, verbose=FALSE) {

  if (missing(pos)) pos = NULL  
  if (!is.null(pos) && !is.numeric(pos)) stop("pos must be numeric.")
  if (!is.null(pos) && is.unsorted(pos)) stop("pos must be in increasing order.")
  if (!is.null(pos) && any(diff(pos)==0)) stop("pos must contain distinct values.")
  if (missing(X)) X = NULL
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  if (!is.null(pos) && is.null(X) && length(pos)!=length(y)) stop("Dimensions don't match [length(pos) != length(y)].")
  if (!is.null(pos) && !is.null(X) && length(pos)!=ncol(X)) stop("Dimensions don't match [length(pos) != ncol(X)].")
  
  if (gamma!=0) {
    nv = if (is.null(X)) length(y) else ncol(X)
    D = getD1dSparse(nv)
    out = fusedlasso(y,X,D,NULL,gamma,approx,maxsteps,minlam,rtol,btol,eps,verbose)
    out$ord = 0
    out$pos = pos
  }
  else out = trendfilter(y,pos,X,0,approx,maxsteps,minlam,rtol,btol,eps,verbose)

  out$call = match.call()
  class(out) = c("fusedlasso", "trendfilter", "genlasso", "list")
  return(out)
}
