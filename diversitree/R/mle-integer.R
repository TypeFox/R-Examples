## Support for optimising one dimensional integer functions.

## All I have implemented here is a simple golden bisection search,
## but there are a number of support functions.  The reason for doing
## bisection rather than one of the quadratic methods is that I am
## concerned that the bounds may not collapse in nicely with quadratic
## search.
do.mle.search.int1d <- function(func, x.init, control, lower, upper) {
  if ( length(x.init) != 1 )
    stop("'int1d' can only be used on univariate problems")
  x.init <- check.integer(x.init)

  func2 <- invert(func)

  control <- modifyList(list(y.init=NULL, interval=NULL, step=1L,
                             factor=2L, maxiter=50, cache=TRUE),
                        control)
  ans <- minimise.int1d(func2, x.init, -control$y.init,
                        control$interval, lower, upper, control)
  ans$value <- -ans$value
  names(ans)[names(ans) == "value"] <- "lnLik"

  ans
}

## Simple wrapper for performing minimisation, rather than
## maximisation - designed to use internally with subplex.mixed()
## (which at this point has a function to minimise).
minimise.int1d <- function(func, x.init, y.init=NULL,
                           interval=NULL, lower=-Inf, upper=Inf,
                           control=list(), ...) {
  if ( length(x.init) != 1 )
    stop("'int1d' can only be used on univariate problems")
  x.init <- check.integer(x.init)
  if ( is.null(y.init) )
    y.init <- func(x.init, ...)
  control <- modifyList(list(step=1, factor=2L, maxiter=50,
                             cache=TRUE), control)

  if ( control$cache )
    func2 <- caching.integer(function(x) func(x, ...), x.init, y.init)
  else
    func2 <- function(x) func(x, ...)

  if ( is.finite(lower) ) lower <- check.integer(lower)
  if ( is.finite(upper) ) upper <- check.integer(upper)
  factor <- check.integer(control$factor)
  step <- check.integer(control$step)
  x.init <- check.integer(x.init)

  if ( is.null(interval) ) {
    tmp <- bracket.int1d(func2, x.init, y.init,
                         control$step, control$factor,
                         control$maxiter, lower, upper)
    interval <- tmp$x[c(1,3)]
    x.init <- tmp$x[2]
    y.init <- tmp$y[2]
  } else {
    interval <- check.integer(sort(interval))
    check.bounds(interval[1], interval[2], x.init)
  }

  ans <- bisect.int1d(func2, x.init, interval)
  ret <- list(par=ans[[1]], value=ans[[2]])
  if ( control$cache )
    ret$count <- length(environment(func2)$cache.x) - 1

  ret
}

## Slightly modified bracketing function, for integer-argument
## functions.  This just adjusts some arguments so that all relevant
## arguments are integers (this prevents ever proposing a non-integer
## argument).
bracket.int1d <- function(func, x.init, y.init, step, factor, maxiter,
                          lower=-Inf, upper=Inf) {
  bracket.1d(func, x.init, y.init, step, factor, maxiter, lower,
             upper)
}

## The functions here are not written particularly nicely, but should
## get the job done.  There is not much effort to reduce the number of
## function evaluations, as it is easiest to use a caching function
## (as caching.integer() makes).
bisect.int1d <- function(func, x.init, interval) {
  R <- 2/(1 + sqrt(5)) # 0.61803399
  C <- 1 - R

  x0 <- interval[1]
  x3 <- interval[2]

  if ( abs(x3-x.init) > abs(x.init-x0) ) { 
    x1 <- x.init
    x2 <- x.init+round(C*(x3-x.init))
  } else { 
    x2 <- x.init
    x1 <- x.init-round(C*(x.init-x0))
  }

  f1 <- func(x1)
  f2 <- func(x2)

  while ( abs(x3-x0) > 1 ) {
    if( f2 < f1 ) {
      x0 <- x1
      x1 <- x2
      x2 <- round(R*x1+C*x3)
      f1 <- f2
      f2 <- func(x2)
    } else {
      x3 <- x2
      x2 <- x1
      x1 <- round(R*x2+C*x0)
      f2 <- f1
      f1 <- func(x1)
    }
  }

  if ( f1 < f2 )
    c(x1, f1)
  else
    c(x2, f2)
}

## This function wraps a function whose first argument is an integer
## (and is the only thing assumed to vary!); if the argument has been
## seen before, the cached value is returned, otherwise the function
## is evaluated and the computed version is both cached and returned.
## This will be suitable for any function where the function
## evaluation is nontrivial.
##
## If initial values are already known (one point, or a series of
## points) these can be passed in as 'x' and 'y', but this is
## optional.
caching.integer <- function(f, x=integer(), y=numeric()) {
  cache.x <- x
  cache.y <- y
  ret <- function(...) {
    i <- as.integer(..1)
    j <- match(i, cache.x)
    if ( is.na(j) ) {
      ans <- f(...)
      cache.x <<- c(cache.x, i)
      cache.y <<- c(cache.y, ans)
    } else {
      ans <- cache.y[[j]]
    }
    ans
  }
}
