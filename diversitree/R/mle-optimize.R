## Line minimisation:
do.mle.search.optimize <- function(func, x.init, control, lower,
                                   upper) {
  func2 <- invert(func)

  if ( length(x.init) != 1 )
    stop("'optimize' can only be used on univariate problems")

  control <- modifyList(list(y.init=NULL, interval=NULL, step=1,
                             maxiter=50, factor=(1 + sqrt(5))/2,
                             useR=FALSE, tol=.Machine$double.eps^0.25),
                        control)

  y.init <- control$y.init
  if ( is.null(y.init) && (control$useR || is.null(control$interval)) )
    y.init <- func2(x.init)
 
  if ( is.null(control$interval) ) {
    tmp <- bracket.1d(func2, x.init, y.init, control$step,
                      control$factor, control$maxiter, lower,
                      upper)
    interval <- tmp$x[c(1,3)]
    x.init <- tmp$x[2]
    y.init <- tmp$y[2]
  } else {
    interval <- control$interval
  }

  if ( control$useR )
    ## Here, we should use the R version, that allows us to recycle
    ## know function evaluations, both in and out of the function.
    ## However, this does not work well.
    stop("Not yet implemented")
  else
    ans <- optimize(func2, interval, maximum=FALSE, tol=control$tol)

  list(par=ans[[1]], lnLik=as.numeric(ans[[2]]))
}

## There is quite a lot of book-keeping here, but this is a simple
## bracketing function, inspired by the one in NR, but with bounds
## added.

## Assume for now that step is positive, so that we have the ordered
## points [x0,x1]
##
##   * If y0 < y1, the optimum lies to the right of y1, so keep the
##     sign of step.
##
##   * If y0 > y1, the optimum lies to the left of y1, so flip the
##     sign of step to search backwards.

## This brackets a *minimum* - beware!
bracket.1d <- function(func, x0, y0, step, factor, maxiter, lower=-Inf,
                    upper=Inf) {
  newp <- function(x, step) {
    xnew <- x + step
    at.edge <- FALSE    
    if ( step > 0 && xnew > upper ) {
      xnew <- upper
      at.edge <- TRUE
    } else if ( step < 0 && xnew < lower ) {
      xnew <- lower
      at.edge <- TRUE
    }
    c(xnew, func(xnew), at.edge)
  }

  xy0 <- c(x0, y0, FALSE)
  xy1 <- newp(x0, step)

  up <- xy0[2] > xy1[2]
  if ( !up ) {
    step <- -step
    tmp <- xy1
    xy1 <- xy0
    xy0 <- tmp
  }
  xy2 <- newp(xy1[1], step * factor)

  for ( i in seq_len(maxiter) ) {
    if ( xy2[3] || xy2[2] > xy1[2] ) {
      m <- rbind(xy0, xy1, xy2)
      return(list(x=m[,1], y=m[,2], at.edge=m[,3]))
    }
    
    xy3 <- newp(xy2[1], step * factor^(i+1))
    xy0 <- xy1
    xy1 <- xy2
    xy2 <- xy3
  }
  stop("Failed to isolate edge - increase maxiter or step?")
}  

## There is a problem with brent() here - in that the initial point
## should not probably be where it is (the termination test may
## trigger at the beginning?).  Consider fixing this?

## brent <- function(f, x.init, interval, y.init, tol) {
##   ## Squared inverse of the golden ratio
##   const <- (3. - sqrt(5.)) * .5

##   tol1 <- 1 + .Machine$double.eps # smallest 1.000... > 1
##   tol3 <- tol / 3.0
##   eps <- .Machine$double.eps^0.25

##   a <- interval[1]
##   b <- interval[2]
##   x <- w <- v <- x.init

##   d <- e <- 0
##   fx <- fw <- fv <- y.init

##   repeat {
##     xm <- (a + b) * 0.5
    
##     tol1 <- eps * abs(x) + tol3
##     t2 <- tol1 * 2.0

##     ## Check stopping criterion
##     if ( abs(x - xm) <= t2 - (b - a) * .5 )
##       break
##     p <- q <- r <- 0.0

##     if ( abs(e) > tol1 ) { # fit parabola
##       r <- (x - w) * (fx - fv)
##       q <- (x - v) * (fx - fw)
##       p <- (x - v) * q - (x - w) * r
##       q <- (q - r) * 2.
##       if (q > 0.) p <- -p else q <- -q
##       r <- e
##       e <- d
##     }

##     if (abs(p) >= abs(q * .5 * r) ||
##         p <= q * (a - x) ||
##         p >= q * (b - x)) { # a golden-section step
##       e <- if (x < xm) b - x else a - x
##       d <- round(const * e)
##     } else { # parabolic-interpolation step
##       d <- p / q
##       u <- x + d

##       if (u - a < t2 || b - u < t2) {
## 	d <- tol1
## 	if (x >= xm)
##           d <- -d
##       }
##     }

##     ## f must not be evaluated too close to x
##     u <- round(x + d)
##     fu <- f(u)

##     ## update  a, b, v, w, and x
##     if (fu <= fx) {
##       if (u < x) b <- x else a <- x
##       v <- w;    w <- x;   x <- u;
##       fv <- fw; fw <- fx; fx <- fu;
##     } else {
##       if (u < x) a <- u else b <- u
##       if (fu <= fw || w == x) {
## 	v <- w; fv <- fw;
## 	w <- u; fw <- fu;
##       } else if (fu <= fv || v == x || v == w) {
## 	v <- u; fv <- fu;
##       }
##     }
##   }

##   ##  end of main loop
##   c(x, fx)
## }

