
## =============================================================================
## hessian    : generates the hessian matrix by numerical differencing
## =============================================================================

hessian <- function (f, x, centered=FALSE, pert=1e-8,  ...)
{
  if (!is.numeric(x))
      stop("x-values should be numeric")
  refx <- x
  reff <- gradient(f,x, ...)
  Nx <- length(x)
  Nf <- length(reff)
  delt <- perturb(x,pert)
  hess <- matrix(nrow = Nf, ncol = Nx, data = 0)
  for (j in 1:Nx) {
    x[j] <- x[j] + delt[j]
    newf <- gradient(f,x, centered=centered, ...)
    del <- (newf - reff)/delt[j]
    if (centered)  {
      # backward formula
      x[j] <- refx[j]-delt[j]
      # recalculate model function value
      newf  <- gradient(f,x, centered=centered, ...)
      del   <- (del-(newf-reff)/delt[j])/2
    }

    hess[, j] <- del
    x[j] <- refx[j]
  }
  colnames(hess) <- names(x)
  rownames(hess) <- attr(del, "names")
  return(hess)
}
