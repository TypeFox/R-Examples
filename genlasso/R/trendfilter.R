# We compute the solution path of the trend filtering problem:
#
# \hat{\beta}(\lambda) =
# \argmin_\beta \|y - X \beta|_2^2 + \lambda\|D \beta\|_1,
#
# where D is (p-k-1) x p is the discrete difference operator of
# order k+1, and X is n x p and full column rank. The solution is
# a piecewise polynomial of degree k, with adaptively chosen knots.

trendfilter <- function(y, pos, X, ord=1, approx=FALSE, maxsteps=2000,
                        minlam=0, rtol=1e-7, btol=1e-7, eps=1e-4, verbose=FALSE) {

  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (missing(pos)) pos = NULL
  if (!is.null(pos) && !is.numeric(pos)) stop("pos must be numeric.")
  if (!is.null(pos) && is.unsorted(pos)) stop("pos must be in increasing order.")
  if (!is.null(pos) && any(diff(pos)==0)) stop("pos must contain distinct values.")
  if (missing(X)) X = NULL
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  if (!is.null(pos) && is.null(X) && length(pos)!=length(y)) stop("Dimensions don't match [length(pos) != length(y)].")
  if (!is.null(pos) && !is.null(X) && length(pos)!=ncol(X)) stop("Dimensions don't match [length(pos) != ncol(X)].")
  if (ord<0 || round(ord)!=ord) stop("ord must be a nonnegative integer.")
  if (length(y) <= ord+1) stop("Not enough data points to fit a trend of the given order [need length(y) > ord+1].")
  if (ord>3) warning(paste("For numerical stability, it is not recommended to run",
                           "trend filtering with a polynomial order larger than 3."))

  # For simplicity
  y = as.numeric(y)
  
  if (is.null(X)) {
    n = length(y)
    if (is.null(pos)) D = getDtfSparse(n,ord)
    else D = getDtfPosSparse(n,ord,pos)
    out = dualpathWideSparse(y,D,NULL,approx,maxsteps,minlam,rtol,btol,verbose)

    # Compute beta, fit, y, bls (this would have been done by
    # the dualpath function)
    out$beta = as.matrix(y - t(D)%*%out$u)
    colnames(out$beta) = colnames(out$u)
    out$fit = out$beta
    out$y = y
    out$bls = y
    
    # Hijack the pathobjs component (to introduce what would
    # have been put here from the dualpath function)
    out$pathobjs$n0 = n
    out$pathobjs$y0 = y
    out$pathobjs$j = 1:n
    out$pathobjs$D0 = D
    out$pathobjs$coldif = 0
  }

   else {
    n = nrow(X)
    p = ncol(X)
    if (length(y)!=n) stop("Dimensions don't match [length(y) != nrow(X)].")

    # Figure out whether we are adding a ridge penalty
    if (p <= n) eps = 0 
    else warning(sprintf(paste("Adding a small ridge penalty (multiplier %g),",
                               "because X has more columns than rows."),eps))
    
    if (is.null(pos)) D = getDtfSparse(p,ord)
    else D = getDtfPosSparse(p,ord,pos)
    out = dualpathTrendX(y,pos,X,D,ord,approx,maxsteps,minlam,rtol,btol,eps,verbose)
  }
  
  out$ord = ord
  out$pos = pos
  out$call = match.call()
  
  class(out) = c("trendfilter", "genlasso", "list")
  return(out)
}
