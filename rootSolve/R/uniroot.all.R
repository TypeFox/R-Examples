
## =============================================================================
## uniroot.all: multiple roots of one nonlinear equation
## =============================================================================

uniroot.all <- function (f, interval, lower= min(interval),
        upper= max(interval), tol= .Machine$double.eps^0.2,
        maxiter= 1000, n = 100, ... ) {

## error checking as in uniroot...
  if (!missing(interval) && length(interval) != 2)
     stop("'interval' must be a vector of length 2")
  if (!is.numeric(lower) || !is.numeric(upper) || lower >=
     upper)
    stop("lower < upper  is not fulfilled")
 
## subdivide interval in n subintervals and estimate the function values
  xseq <- seq(lower,upper,len=n+1)
  mod  <- f(xseq,...)

## some function values may already be 0
  Equi <- xseq[which(mod==0)]

  ss   <- mod[1:n]*mod[2:(n+1)]  # interval where functionvalues change sign
  ii   <- which(ss<0)

  for (i in ii)
    Equi <- c(Equi,uniroot(f,lower=xseq[i],upper=xseq[i+1],...)$root)

  return(Equi)
}

