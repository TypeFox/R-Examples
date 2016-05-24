# We compute the solution path of the generalized lasso problem:
#
# \hat{\beta}(\lambda) =
# \argmin_\beta \|y - X \beta|_2^2 + \lambda\|D \beta\|_1,
#
# where X is n x p and D is m x p. The penalty matrix D can
# be arbitrary, the predictor matrix X should have rank(X)=p;
# if not, we regularize by adding a small ridge penalty.

genlasso <- function(y, X, D, approx=FALSE, maxsteps=2000, minlam=0,
                     rtol=1e-7, btol=1e-7, eps=1e-4, verbose=FALSE,
                     svd=FALSE) {
 
  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 1].")
  if (missing(X)) X = NULL
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  if (missing(D)) stop("D is missing.")
  if (!is.matrix(D) && c(attributes(class(D))$package,"")[[1]] != "Matrix") {
    stop("D must be a matrix or a Matrix (from the Matrix package).")
  }
  if (is.null(X) && length(y)!=ncol(D)) stop("Dimensions don't match [length(y) != ncol(D)].")
  if (checkrows(D)) stop("D cannot have duplicate rows.")

  # For simplicity
  y = as.numeric(y)
  
  # X should be treated as the identity
  if (is.null(X)) {
    if (!svd) out = dualpath(y,D,approx,maxsteps,minlam,rtol,btol,verbose)
    else out = dualpathSvd(y,D,approx,maxsteps,minlam,rtol,btol,verbose)
  }
  
  # X is given
  else {
    if (!is.matrix(D)) {
      warning("Converting D to a dense matrix, because X is not the identity.")
      D = as.matrix(D)
    }

    n = nrow(X)
    p = ncol(X)
    if (length(y)!=n) stop("Dimensions don't match [length(y) != nrow(X)].")
    if (ncol(D)!=p) stop("Dimensions don't match [ncol(X) != ncol(D)].")

    ridge = FALSE
    if (p > n) {
      if (eps<=0) stop("eps must be positive when X has more columns than rows.")
      warning(sprintf("Adding a small ridge penalty (multiplier %g), because X has more columns than rows.",eps))
      ridge = TRUE
    }
    else {
      # Check that X has full column rank
      x = svd(X)
      if (all(x$d >= rtol)) {
        y2 = as.numeric(x$u %*% t(x$u) %*% y)
        Xi = x$v %*% (t(x$u) / x$d)
        D2 = D %*% Xi
      }
      else {
        if (eps<=0) stop("eps must be positive when X is column rank deficient.")
        warning(sprintf("Adding a small ridge penalty (multiplier %g), because X is column rank deficient.",eps))
        ridge = TRUE
      }
    }
    
    if (ridge) {
      x = svd(rbind(X,diag(sqrt(eps),p)))
      y2 = as.numeric(x$u %*% t(x$u) %*% c(y,rep(0,p)))
      Xi = x$v %*% (t(x$u) / x$d)
      D2 = D %*% Xi
    }

    if (!svd) out = dualpath(y2,D2,approx,maxsteps,minlam,rtol,btol,verbose)
    else out = dualpathSvd(y2,D2,approx,maxsteps,minlam,rtol,btol,verbose)

    # Save these path objects for internal use later
    out$pathobjs$y2 = y2
    out$pathobjs$Xi = Xi

    # Fix beta, fit, y, bls, and save the X matrix
    out$beta = Xi %*% out$fit
    out$fit = X %*% out$beta
    out$y = y
    out$bls = Xi %*% y2
    out$X = X
  }

  out$call = match.call()
  class(out) = c("genlasso", "list")
  return(out)
}
