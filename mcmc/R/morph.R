euclid.norm <- function(x) {
  sqrt(sum(x * x))
}

isotropic <- function(f) {
  force(f)
  function(x) {
    x.norm <- euclid.norm(x)
    if (x.norm == 0)
      rep(0, length(x))
    else
      f(x.norm) * x / x.norm
  }
}

isotropic.logjacobian <- function(f, d.f) {
  force(f)
  force(d.f)
  function(x) {
    x.norm <- euclid.norm(x)
    k <- length(x)
    if (x.norm == 0) {
      k * log(d.f(x.norm))
    } else {
      log(d.f(x.norm)) + (k - 1) * (log(f(x.norm)) - log(x.norm))
    }
  }
}

# rearranged by charlie
# add starting point so we can use a better starting point for
#     the current actual use
# also do one more Newton step after error is < sqrt(machine.eps)
#     so ultimate error is about machine.eps
newton.raphson <- function(f, df, x, r, x0 = 2 * r) {
  f.err <- function(cur) f(cur) - x
  step <- function(cur, err=NULL)
    cur - ifelse(is.null(err), f.err(cur), err) / df(cur)
  err <- f.err(x0)
  while(err >= sqrt(.Machine$double.eps)) {
    x0 <- step(x0, err)
    err <- f.err(x0)
  }
  # if you don't want to use an extra-step, replace this return with
  # return(x0)
  return(step(x0, err))
}

subexponential <- function(b=1) {
  if (missing(b) | is.null(b)) b <- 1
  stopifnot(b > 0)
  force(b)
  f.inv <- function(x) ifelse(x > 1/b,
                              exp(b * x) - exp(1)/3,
                              (x * b)^3 * exp(1) / 6 + x * b * exp(1) / 2)
  d.f.inv <- function(x) ifelse(x > 1/b,
                                b * exp(b * x),
                                b * (x * b)^2 * exp(1) / 2 + b * exp(1) / 2)
  f <- function(x) {
    # x > exp(b * 1 / b) - exp(1) / 3
    if (x > 2 * exp(1) / 3) {
      log(x + exp(1)/3) / b
    } else {
      poly.inv <- exp(1/3) *
        (sqrt(b^12 * (9 * x^2 + exp(2))) - 3 * b^6 * x)^(-1/3)
      poly.inv * b - 1 / (poly.inv * b^3)
    }
  }
  return(list(f=f, f.inv=f.inv, d.f.inv=d.f.inv))
}

exponential <- function(r=1, p=3) {
  if (missing(p) || is.null(p)) p <- 3
  if (missing(r) || is.null(r)) r <- 0
  stopifnot(p > 2)
  stopifnot(r >= 0)
  f.inv <- function(x) ifelse(x <= r, x, x + (x-r)^p)
  d.f.inv <- function(x) ifelse(x <= r, 1, 1 + p * (x-r)^(p-1))
  if (p == 3) {
    g <- function(x) {
      n <- sqrt((27*r-27*x)^2 + 108) + 27 * (r - x)
      r + (2/n)^(1/3) - (n/2)^(1/3)/3
    }
    f <- function(x) ifelse(x < r, x, g(x))
  } else {
    # No general closed form solution exists.  However, since the
    # transformation has polynomial form, using the Newton-Raphson method
    # should work well.
    f <- function(x) ifelse(x < r, x,
        newton.raphson(f.inv, d.f.inv, x, r, x0 = r + x^(1 / p)))
  }
  return(list(f=f, f.inv=f.inv, d.f.inv=d.f.inv))
}

.make.outfun <- function(out) {
  force(out)
  function(f) {
    force(f)
    if (is.null(f))
      return(out$inverse)
    else if (is.function(f))
      return(function(state, ...) f(out$inverse(state), ...))
    else
      return(function(state) out$inverse(state)[f])
  }
}

identity.func <- function(x) x
morph.identity <- function() {
  out <- list(transform=identity.func,
              inverse=identity.func,
              lud=function(f) function(x, ...) f(x, ...),
              log.jacobian=function(x) 0,
              center=0,
              f=identity.func,
              f.inv=identity.func)
  out$outfun <- .make.outfun(out)
  return(out)
}

morph <- function(b, r, p, center) {
  if (all(missing(b), missing(r), missing(p), missing(center)))
    return(morph.identity())
  if (missing(center)) center <- 0
  use.subexpo <- !missing(b)
  use.expo <- !(missing(r) && missing(p))
  
  if (!use.expo && !use.subexpo) {
    f <- function(x) x
    f.inv <- function(x) x
    log.jacobian <- function(x) 0
  } else {
    if (use.expo && !use.subexpo) {
      expo <- exponential(r, p)
      
      f <- expo$f
      f.inv <- expo$f.inv
      d.f.inv <- expo$d.f.inv
    } else if (!use.expo && use.subexpo) {
      subexpo <- subexponential(b)
      
      f <- subexpo$f
      f.inv <- subexpo$f.inv
      d.f.inv <- subexpo$d.f.inv
    } else { #use.expo && use.subexpo
      expo <- exponential(r, p)
      subexpo <- subexponential(b)
      
      f <- function(x) expo$f(subexpo$f(x))
      f.inv <- function(x) subexpo$f.inv(expo$f.inv(x))
      d.f.inv <- function(x) expo$d.f.inv(x) * subexpo$d.f.inv(expo$f.inv(x))
    }
    
    f <- isotropic(f)
    f.inv <- isotropic(f.inv)
    log.jacobian <- isotropic.logjacobian(f.inv, d.f.inv)
  }

  out <- list(f=f, f.inv=f.inv, log.jacobian=log.jacobian,
              center=center)
  out$transform <- function(state) out$f(state - out$center)
  out$inverse <- function(state) out$f.inv(state) + out$center
  
  out$outfun <- .make.outfun(out)

  out$lud <- function(lud) {
    force(lud)
    function(state, ...) {
      foo <- lud(out$inverse(state), ...)
      if (length(foo) != 1)
          stop("log unnormalized density function returned vector not scalar")
      if (is.na(foo))
          stop("log unnormalized density function returned NA or NaN")
      if (foo == -Inf) return(foo)
      if (! is.finite(foo))
          stop("log unnormalized density function returned +Inf")
      foo + out$log.jacobian(state)
    }
  }

  return(out)
}

