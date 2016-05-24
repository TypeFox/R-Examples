#' @name internal
#' 
#' @inheritParams kden
#' @inheritParams bckden
#' @inheritParams itmweibullgpd
#' @inheritParams itmgng
#' @inheritParams psden
#' 
#' @title Internal Functions
#'
#' @description Internal functions not designed to be used directly, but are all exported
#' to make them visible to users.
#'
#' @details Internal functions not designed to be used directly. No error
#' checking of the inputs is carried out, so user must be know what they are doing.
#' They are undocumented, but are made visible to the user.
#' 
#' Mostly, these are used in the kernel density estimation functions.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @section Acknowledgments: Based on code
#' by Anna MacDonald produced for MATLAB.
#' 
#' @seealso \code{\link[stats:density]{density}}, \code{\link[evmix:kden]{kden}}
#' and \code{\link[evmix:bckden]{bckden}}.
#' 
NULL

# The following functions are supposed to be used internally for functions
# and are not designed to be used directly by users, there is no error checking


#' @export
#' @rdname internal
kdenx <- function(x, kerncentres, lambda, kernel = "gaussian") {
  mean(kdz((x - kerncentres)/lambda, kernel))/lambda
}

#' @export
#' @rdname internal
pkdenx <- function(x, kerncentres, lambda, kernel = "gaussian") {
  mean(kpz((x - kerncentres)/lambda, kernel))
}

#' @export
#' @rdname internal
bckdenxsimple <- function(x, kerncentres, lambda, kernel = "gaussian") {
  # use notation of Jones (1993)

  # distance of evaluation point to lower boundary (p in paper)
  truncpoint = x/lambda
  
  # distance of evaluation point to kernel centres
  u = (x - kerncentres)/lambda

  # apply adjustment only if point is close to zero boundary
  maxp = ifelse(kernel == "gaussian", 5, 1)
  if (truncpoint <= maxp) {
    # integral of moments upto p
    a0 = ka0(truncpoint, kernel = kernel)
    a1 = ka1(truncpoint, kernel = kernel) 
    a2 = ka2(truncpoint, kernel = kernel) 
  
    # weights in local linear fitting
    denom = (a2*a0 - a1^2)
    lx = a2/denom
    mx = a1/denom
    d = mean((lx - mx*u)*kdz(u, kernel = kernel))/lambda    
  } else {
    d = mean(kdz(u, kernel = kernel))/lambda  
  }

  d
}

#' @export
#' @rdname internal
pbckdenxsimple <- function(x, kerncentres, lambda, kernel = "gaussian") {
    
  # apply adjustment only if point is close to zero boundary
  maxp = ifelse(kernel == "gaussian", 5, 1)*lambda
  if (x <= maxp) {
    # Reuse density function to do numerical integration
    bckdenint = try(integrate(Vectorize(bckdenxsimple, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, kernel = kernel,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of simple boundary corrected KDE")
    }
    p = bckdenint$value    
  } else {
    # Reuse density function to do numerical integration upto 5*lambda
    bckdenint = try(integrate(Vectorize(bckdenxsimple, "x"), lower = 0, upper = maxp,
      kerncentres = kerncentres, lambda = lambda, kernel = kernel,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of simple boundary corrected KDE")
    }
    p = bckdenint$value    

    # Now add contribution from maxp upwards
    u = (x - kerncentres)/lambda
    
    p = p + mean(kpz(u, kernel = kernel) - kpz((maxp - kerncentres)/lambda, kernel = kernel))   
  }
  
  p
}

#' @export
#' @rdname internal
bckdenxcutnorm <- function(x, kerncentres, lambda, kernel = "gaussian") {
  
  # distance of kernel centres to lower boundary
  truncpoint = kerncentres/lambda
  
  # how much of kernel is in range of support, so (1-a0) gives leakage past boundary
  a0 = kpz(truncpoint, kernel)
  
  u = (x - kerncentres)/lambda
  
  mean(kdz(u, kernel)/a0)/lambda
}

#' @export
#' @rdname internal
pbckdenxcutnorm <- function(x, kerncentres, lambda, kernel = "gaussian") {
  
  # distance of kernel centres to lower boundary
  truncpoint = kerncentres/lambda
  
  # how much of kernel is in range of support, so (1-a0) gives leakage past boundary
  a0 = kpz(truncpoint, kernel)
  
  u = (x - kerncentres)/lambda
  
  mean((kpz(u, kernel) - kpz(-kerncentres/lambda, kernel))/a0)
}

#' @export
#' @rdname internal
bckdenxrenorm <- function(x, kerncentres, lambda, kernel = "gaussian") {
  
  # distance of x to lower boundary
  truncpoint = x/lambda
  
  u = (x - kerncentres)/lambda

  maxp = ifelse(kernel == "gaussian", 5, 1)
  if (truncpoint < maxp) {
    # first order correction
    a0 = kpz(truncpoint, kernel)
    
    d = mean(kdz(u, kernel)/a0)/lambda
  } else {
    d = mean(kdz(u, kernel))/lambda    
  }
  d
}

#' @export
#' @rdname internal
pbckdenxrenorm <- function(x, kerncentres, lambda, kernel = "gaussian") {
  
  maxp = ifelse(kernel == "gaussian", 5, 1)*lambda
  if (x <= maxp) {
    # Reuse density function to do numerical integration
    bckdenint = try(integrate(Vectorize(bckdenxrenorm, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, kernel = kernel,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of renormalisation boundary corrected KDE")
    }
    p = bckdenint$value    
  } else {
    # Reuse density function to do numerical integration upto 5*lambda
    bckdenint = try(integrate(Vectorize(bckdenxrenorm, "x"), lower = 0, upper = maxp,
      kerncentres = kerncentres, lambda = lambda, kernel = kernel,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of renormalisation boundary corrected KDE")
    }
    
    p = bckdenint$value

    # Now add contribution from maxp upwards
    u = (x - kerncentres)/lambda
  
    p = p + mean(kpz(u, kernel = kernel) - kpz((maxp - kerncentres)/lambda, kernel = kernel))
  }
  p
}

#' @export
#' @rdname internal
bckdenxreflect <- function(x, kerncentres, lambda, kernel = "gaussian") {
  mean(kdz((x - kerncentres)/lambda, kernel = kernel) + kdz((x + kerncentres)/lambda, kernel = kernel))/lambda
}

#' @export
#' @rdname internal
pbckdenxreflect <- function(x, kerncentres, lambda, kernel = "gaussian") {
  mean(kpz((x - kerncentres)/lambda, kernel = kernel) - kpz(-kerncentres/lambda, kernel = kernel) 
    + kpz(-kerncentres/lambda, kernel = kernel) - kpz(-(x+kerncentres)/lambda, kernel = kernel) )
}

#' @export
#' @rdname internal
pxb <- function(x, lambda) 2*lambda^2 + 2.5 - sqrt(4*lambda^4 + 6*lambda^2 + 2.25 - x^2 - x/lambda)

#' @export
#' @rdname internal
bckdenxbeta1 <- function(x, kerncentres, lambda, xmax) {
  
  # rescale x and kernel centres, then treat as bounded on [0, 1] so beta kernel used directly
  x = x/xmax
  kerncentres = kerncentres/xmax
  lambda = lambda/xmax
  
  if ((x >= 0) & (x <= 1)){
    d = mean(dbeta(kerncentres, x/lambda + 1, (1 - x)/lambda + 1))
  } else {
    d = 0
  }
  d/xmax
}

#' @export
#' @rdname internal
pbckdenxbeta1 <- function(x, kerncentres, lambda, xmax) {
  if (x <= 0) {
    p = 0
  } else if (x > xmax) {
    p = 1
  } else {
    # Re-use density function to do numerical integration
    bckdenint = try(integrate(Vectorize(bckdenxbeta1, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, xmax = xmax,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of beta based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
bckdenxbeta2 <- function(x, kerncentres, lambda, xmax) {
  
  # rescale x and kernel centres, then treat as bounded on [0, 1] so beta kernel used directly
  x = x/xmax
  kerncentres = kerncentres/xmax
  lambda = lambda/xmax

  # split x into cases: middle, near 0, near 1 or outside range
  if ((x >= 2*lambda) & (x <= (1 - 2*lambda))) {
    d = mean(dbeta(kerncentres, x/lambda, (1 - x)/lambda))
  } else if ((x >= 0) & (x < 2*lambda)) {
    d = mean(dbeta(kerncentres, pxb(x, lambda), (1 - x)/lambda))
  } else if ((x > (1 - 2*lambda)) & (x <= 1)) {
    d = mean(dbeta(kerncentres, x/lambda, pxb(1 - x, lambda)))
  } else {
    d = 0
  }
  d/xmax
}

#' @export
#' @rdname internal
pbckdenxbeta2 <- function(x, kerncentres, lambda, xmax) {
  if (x <= 0) {
    p = 0
  } else if (x > xmax) {
    p = 1
  } else {
    # Reuse density function to do numerical integration
    bckdenint = try(integrate(Vectorize(bckdenxbeta2, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, xmax = xmax,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of beta based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
bckdenxgamma1 <- function(x, kerncentres, lambda) {
  
  if (x >= 0) {
    d = mean(dgamma(kerncentres, shape = x/lambda + 1, scale = lambda))
  } else {
    d = 0
  }
  d
}

#' @export
#' @rdname internal
pbckdenxgamma1 <- function(x, kerncentres, lambda) {

  if (x < 0) {
    p = 0
  } else {
    # Reuse density function to do numerical integration
    bckdenint = try(integrate(Vectorize(bckdenxgamma1, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda,
      subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
  
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of gamma based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
bckdenxgamma2 <- function(x, kerncentres, lambda) {

  if ((x >= 0) & (x < 2*lambda)) {
    d = mean(dgamma(kerncentres, shape = (x/lambda)^2/4 + 1, scale = lambda))
  } else if (x >= 2*lambda) {
    d = mean(dgamma(kerncentres, shape = x/lambda, scale = lambda))
  } else {
    d = 0
  }
  d
}

#' @export
#' @rdname internal
pbckdenxgamma2 <- function(x, kerncentres, lambda) {

  if (x < 0) {
    p = 0
  } else {
    # Reuse density function to do numerical integration
    bckdenint = try(integrate(Vectorize(bckdenxgamma2, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda,
      subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
  
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of gamma based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
bckdenxcopula <- function(x, kerncentres, lambda, xmax) {
  # rescale x and kernel centres, then treat as bounded on [0, 1] so beta kernel used directly
  x = x/xmax
  kerncentres = kerncentres/xmax
  
  # quantities used in Gaussian copula kernels
  lambda2 = lambda^2
  stdinv = ifelse((x >= 0) & (x <= 1), qnorm(x), NA)
  kerninv = ifelse((kerncentres >= 0) & (kerncentres <= 1), qnorm(kerncentres), NA)
  
  # Renormalisation constant and (unscaled) density
  d1 = exp(-((1 - lambda2)*stdinv)^2/2/lambda2/(2 - lambda2))/lambda/sqrt(2 - lambda2)
  d2 = mean(exp(-(1 - lambda2)*((1 - lambda2)*kerninv^2 - 2*stdinv*kerninv)/
      2/lambda2/(2 - lambda2)))
  
  # treat those outside range [0, xmax] as 0 density
  d = ifelse((x >= 0) & (x <= 1), ifelse(is.infinite(d2), 0, d1*d2/xmax), 0)
}

#' @export
#' @rdname internal
pbckdenxcopula <- function(x, kerncentres, lambda, xmax) {
  if (x <= 0) {
    p = 0
  } else if (x > xmax) {
    p = 1
  } else {
    # Reuse density function to do numerical integration
    bckdenint = try(integrate(Vectorize(bckdenxcopula, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, xmax = xmax,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of copula based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
pbckdenxlog <- function(x, kerncentres, lambda, offset, kernel = "gaussian") {

  # Reuse density function to do numerical integration
  bckdenint = try(integrate(dbckden, lower = 0, upper = x,
    kerncentres = kerncentres, lambda = lambda, kernel = kernel, bcmethod = "logtrans",
    offset = offset, subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
  
  if (inherits(bckdenint, "try-error")) {
    bckdenint$value = NA
    warning("failed to numerically evaluate cdf of log transform based boundary corrected KDE")
  }
  bckdenint$value
}

#' @export
#' @rdname internal
pbckdenxnn <- function(x, kerncentres, lambda, kernel = "gaussian", nn) {

  maxp = ifelse(kernel == "gaussian", 5, 1)*lambda
  # apply adjustment only if point is close to zero boundary
  if (x <= maxp) {
    # Reuse density function to do numerical integration
    bckdenint = try(integrate(dbckden, lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, kernel = kernel,
      bcmethod = "simple", proper = FALSE, nn = nn,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
  
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of non-negative (simple) boundary corrected KDE")
    }
    p = bckdenint$value
  } else {
    # Reuse density function to do numerical integration upto maxp
    bckdenint = try(integrate(dbckden, lower = 0, upper = maxp,
      kerncentres = kerncentres, lambda = lambda, kernel = kernel,
      bcmethod = "simple", proper = FALSE, nn = nn,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
  
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of non-negative (simple) boundary corrected KDE")
    }
    p = bckdenint$value
    
    # Now add contribution from maxp upwards
    u = (x - kerncentres)/lambda
    
    p = p + mean(kpz(u, kernel = kernel) - kpz((maxp - kerncentres)/lambda, kernel = kernel))   
  }
  p
}

#' @export
#' @rdname internal
qmix <- function(x, u, epsilon) {
  # Holden and Haug's q mixing function
  # (negate x and u to get p mixing function)
  
  xu = abs(x - u)
  whichi = which(xu < epsilon)
  whichu = which(x >= (u + epsilon))
  
  qm = x
  qm[whichi] = (x[whichi] + u - epsilon)/2 + epsilon * cos(pi * (x[whichi] - u)/epsilon/2)/pi
  qm[whichu] = u # different, but more computationally convenient, definition (was x-epsilon)
  qm
}

#' @export
#' @rdname internal
qmixprime <- function(x, u, epsilon) {
  # Holden and Haug's q' mixing function
  # (negate x and u to get p' mixing function)
  
  xu = abs(x - u)
  whichi = which(xu < epsilon)
  whichu = which(x >= (u + epsilon))
  
  qprime = rep(1, length(x))
  qprime[whichi] = 1/2 - sin(pi * (x[whichi] - u)/epsilon/2)/2
  qprime[whichu] = 0  # different, but more computationally convenient, definition (was 1)
  qprime
}

#' @export
#' @rdname internal
qgbgmix <- function(x, ul, ur, epsilon) {
  # Holden and Haug's q_i mixing function reformulated for
  # middle bulk model, with GPD for both tails
  
  # Note that intervals must not overlap (not checked)
  midpoint = (ul + ur)/2
  whichl = which(x < midpoint)
  whichu = which(x >= midpoint)
  
  qm = x
  qm[whichl] = -qmix(-x[whichl], -ul, epsilon)
  qm[whichu] = qmix(x[whichu], ur, epsilon)
  qm
}

#' @export
#' @rdname internal
qgbgmixprime <- function(x, ul, ur, epsilon) {
  # Holden and Haug's q'_i mixing function reformulated for
  # middle bulk model, with GPD for both tails
  
  # Note that intervals must not overlap (not checked)  
  midpoint = (ul + ur)/2
  whichl = which(x < midpoint)
  whichu = which(x >= midpoint)
  
  qmprime = x
  qmprime[whichl] = qmixprime(-x[whichl], -ul, epsilon)
  qmprime[whichu] = qmixprime(x[whichu], ur, epsilon)
  qmprime
}

#' @export
#' @rdname internal
pscounts <- function(x, beta, design.knots, degree) {
  bsplines = splineDesign(design.knots, x, ord = degree + 1)
  d = exp(bsplines %*% beta)
  d
}

