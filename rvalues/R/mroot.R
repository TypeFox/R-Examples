mroot <- function(f, lower, upper, ..., f.lower = f(lower, ...), f.upper = f(upper, ...),
                  tol = .Machine$double.eps^0.25, maxiter = 5000) {
  if (!is.numeric(lower) || !is.numeric(upper) || any(lower >= 
        upper))  {
    stop("lower < upper  is not fulfilled")
  }
  if (any(is.na(f.lower)))  {
    stop("f.lower = f(lower) is NA")
  }
  if (any(is.na(f.upper)))  {
    stop("f.upper = f(upper) is NA")
  }
  if(length(lower)!=length(f.lower) | length(upper)!=length(f.upper)) {
     stop("f.lower and lower must have the same length and f.upper
          and upper must have the same length")  
  }
  if (any(f.lower * f.upper > 0)) {
    stop("f() values at end points not of opposite sign")
  }
  if(tol <= 0.0) {
    stop("non-positive tolerance value")
  }
  f.check <- function(x, ...) {
    ff <- f(x, ...)
    as.double(ff)
  }
  st <- length(lower)
  #mid <- 0.5*(lower + upper)
  res <- .Call("ZeroIn", as.double(lower), as.double(upper), as.double(f.lower), as.double(f.upper),
               body(f.check), as.integer(st), as.double(tol), as.integer(maxiter), 
               new.env(), PACKAGE="rvalues");

  nn <- length(lower) + 1
  iter <- res[nn]
  if (iter < 0) {
    warning(sprintf(ngettext(maxiter, "_NOT_ converged in %d iteration", 
                             "_NOT_ converged in %d iterations"), maxiter), domain = NA)
    iter <- maxiter
  }
  ans <- list()
  ans$root <- res[-nn]
  ans$iter <- iter
  return(ans)
}
