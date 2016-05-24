`bvls` <- function(A, b, bl, bu, key = 0, istate=rep(0, ncol(A)+1)) {
 
  M <- nrow(A) 
  N <- ncol(A)
  X <- rep(0, N)

  if(! is.numeric(A)) {
      cat(" 'A' must be numeric \n")
      return() 
  }
  if(! is.numeric(b)) {
      cat(" 'b' must be numeric \n")
      return() 
  }
  if(length(bu)!=N||length(bl)!=N) {
      cat("bounds not correct length \n")
      return() 
  }
  if(! (is.numeric(bu) && is.numeric(bl))) {
      cat("bounds contain non-numeric values \n")
      return() 
  }

  ## working arrays 
  W <- rep(0, N)
  mm <- min(M,N)
  act <- rep(0, M*(mm+2))
  zz  <- rep(0, M)
  loopA <- 0

  sol <- .Fortran("bvls", key = key, m = M, n = N, 
  a = A, b = b, bl = bl, bu = bu, x = X, 
  w=W, act = act, zz = zz, istate = istate, 
  loopA = loopA, PACKAGE="bvls", NAOK=TRUE)

  fitted <- A %*% sol$x 
  resid <- b - fitted 
  bvls.out <- list(x=sol$x, deviance=sum(resid^2),
                   residuals=resid, fitted = fitted)
  class(bvls.out) <- "bvls"
  bvls.out 

}
print.bvls <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("Bounded-variable least squares model\n")

    cat("x estimates:", x$x, "\n")
    cat("residual sum-of-squares: ", format(x$deviance, digits = digits),
	"\n", sep = '')
    
    invisible(x)
}
residuals.bvls <- function(object,...)  object$residuals
coef.bvls <- function(object,...) object$x
fitted.bvls <- function(object,...) object$fitted 
deviance.bvls <- function(object,...) object$deviance

