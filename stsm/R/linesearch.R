
step.maxsize <- function(x, xlo, xup, pd, cap = 1)
{
  # choose lambda such that 
  # x_i + lambda * pd_i in [xlo_i, xup_i]
  # for each parameters 'x_i'

  if (all(pd == 0))
    return(cap)

  ref.lo <- ref.up <- cap

  #ref <- which(pd < 0)
  ref <- intersect(which(xlo != -Inf), which(pd < 0))
  if (length(ref) > 0)
  {
    ref.lo <- min((xlo[ref] - x[ref]) / pd[ref])
    ref.lo <- subset(ref.lo, ref.lo > 0)
  }

  #ref <- which(pd > 0)
  ref <- intersect(which(xup != Inf), which(pd > 0))
  if (length(ref) > 0)
  {
    ref.up <- min((xup[ref] - x[ref]) / pd[ref])
    ref.up <- subset(ref.up, ref.up > 0)
  }

  a <- min(c(ref.lo, ref.up, cap))

#debug
stopifnot(a >= 0)

  a
}

Brent.fmin <- function(a = 0, b, fcn, tol = .Machine$double.eps^0.25, ...)
#
# ported from R sources, procedure Brent_fmin in file fmin.c
#
{
  counts <- c("fcn" = 0, "grd" = NA)

  # squared inverse of the golden ratio
  c <- (3 - sqrt(5)) * 0.5;

  eps <- .Machine$double.eps
  tol1 <- eps + 1 # the smallest 1.000... > 1
  eps <- sqrt(eps)
 
  v <- a + c * (b - a)
  vx <- x <- w <- v

  d <- e <- 0

  fx <- fcn(x, ...)
  counts[1] <- counts[1] + 1

  fw <- fv <- fx

  tol3 <- tol / 3

  # main loop

  iter <- 0
  cond <- TRUE
  while (cond)
  {
    xm <- (a + b) * 0.5
    tol1 <- eps * abs(x) + tol3
    t2 <- tol1 * 2

    # check stopping criterion

    if (abs(x - xm) <= t2 - (b - a) * 0.5)
      break

    r <- q <- p <- 0

    # fit parabola

    if (abs(e) > tol1)
    {
      r <- (x - w) * (fx - fv)
      q <- (x - v) * (fx - fw)
      p <- (x - v) * q - (x - w) * r
      q <- (q - r) * 2
      if (q > 0) {
        p <- -p } else 
        q <- -q
      r <- e
      e <- d
    }

    if (abs(p) >= abs(q * 0.5 * r) ||
      p <= q * (a - x) || p >= q * (b - x))
    { # golden-section step
      if (x < xm) {
        e <- b - x } else 
        e <- a - x

      d <- c * e
    } else { # parabolic-interpolation step 
      d <- p / q
      u <- x + d

      # f must not be evaluated too close to a or b
      if (u - a < t2 || b - u < t2)
      {
        d <- tol1
        if (x >= xm)
          d <- -d
      }
    }

    # f must not be evaluated too close to x

    if (abs(d) >= tol1) {
      u <- x + d } else if (d > 0) {
      u <- x + tol1 } else
      u <- x - tol1

    fu <- fcn(u, ...)
    counts[1] <- counts[1] + 1

    # update  a, b, v, w, and x

    if (fu <= fx)
    {
      if (u < x) { 
        b <- x } else 
        a <- x
      v <- w; w <- x; x <- u
      vx <- c(vx, x)
      fv <- fw; fw <- fx; fx <- fu
    } else {
      if (u < x) {
        a <- u } else 
        b <- u

      if (fu <= fw || w == x)
      {
        v <- w; fv <- fw
        w <- u; fw <- fu
      } else if (fu <= fv || v == x || v == w)
      {
        v <- u
        fv <- fu
      }
    }

    iter <- iter + 1
  }

#return same output in 'minimum 'and 'x', 
#the former to follow same naming as in 'optimize'
  list(vx = vx, minimum = x, x = x, fx = fx, iter = iter, counts = counts)
}

linesearch <- function(b, fcn, grd, ftol = 0.0001, gtol = 0.9, ...)
#
# Based on Nocedal-Wright 2006 chapter 3
# A preliminary version was based on description given in Brent (YYYY)
# may be worth another look at it and the comments
#
{
  approx.quad2 <- function(x0, x1, f0, f1, g0)
  #code taken from linesearch code in /explore/optimization
  {
    dx <- x0 - x1
    dx2 <- x0^2 - x1^2

    a <- (f0 - f1 - g0*dx) / (dx2 - 2*x0*dx)
    b <- g0 - 2*a*x0

    x <- -b / (2*a)

    if (is.na(x) || (x < min(x0, x1)) || (x > max(x0, x1)))
      return ((x0 + x1) / 2)
    x
  }

  wolfe.cond1 <- function(a, fa, f0, g0, ftol)
  {
    fa > f0 + ftol * a * g0
  }

  wolfe.cond2 <- function(ga, mc2.g0)
  {
    abs(ga) <= mc2.g0
  }

  counts = c("fcn" = 0, "grd" = 0)
  cond <- TRUE

  a <- 0
  alpha <- c(a, (b + a) / 2)

  f0 <- fcn(0, ...)
  g0 <- grd(0, ...)
  counts <- counts + c(1, 1)

  mc2.g0 <- -gtol * g0

  # bracketing

  while (cond)
  {
    i <- length(alpha)

    if (i > 2) {
      alpha.old <- ai
      falpha.old <- falpha
      galpha.old <- galpha
    } else {
      alpha.old <- 0
      falpha.old <- f0
      galpha.old <- g0
    }

    ai <- alpha[i]
    falpha <- fcn(ai, ...)
    galpha <- grd(ai, ...)
    counts <- counts + c(1, 1)

    if (wolfe.cond1(ai, falpha, f0, g0, ftol) || 
      (i > 2 && falpha >= falpha.old))
    {
      x <- sort(c(alpha.old, ai), index.return = TRUE)
      y <- c(falpha.old, falpha)[x$ix]
      z <- c(galpha.old, galpha)[x$ix]
      x <- x$x

      break
    }

    if (wolfe.cond2(galpha, mc2.g0))
    {
      cond <- FALSE
      break      
    }

    if (galpha >= 0)
    {
      x <- sort(c(alpha.old, ai), index.return = TRUE)
      y <- c(falpha.old, falpha)[x$ix]
      z <- c(galpha.old, galpha)[x$ix]
      x <- x$x

      break
    }

    tmp <- (ai + b) / 2
    alpha <- c(alpha, tmp)

    if (length(alpha) > 30)
      cond <- FALSE
  }

  # sectioning

  while (cond)
  {
    tmp <- approx.quad2(x[1], x[2], y[1], y[2], z[1])

    alpha <- c(alpha, tmp)

    i <- length(alpha)

    ai <- alpha[i] # tmp
    falpha <- fcn(ai, ...)
    counts <- counts + c(1, 0)

    if (wolfe.cond1(ai, falpha, f0, g0, ftol) || 
      falpha >= y[1]) {
      x[2] <- ai
      y[2] <- falpha
    } else {

      galpha <- grd(ai, ...)
      counts <- counts + c(0, 1)

      if (wolfe.cond2(galpha, mc2.g0))
        break

      if (galpha * (x[2] - x[1]) >= 0)
      {
        x[2] <- x[1]
        y[2] <- y[1]        
      }

      x[1] <- ai
      y[1] <- y[1]
      z[1] <- galpha      
    }

    if (length(alpha) > 50)
      cond <- FALSE
  }

  n <- length(alpha)
  x <- alpha[n]
  list(vx = alpha, minimum = x, x = x, fx = falpha,
    iter = n - 1, counts = counts) #convergence = cond
}
