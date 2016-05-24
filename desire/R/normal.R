##
## normal.R - Desirability functions based on the Normal Distribution
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.tu-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.tu-dortmund.de>
##

##
## normOpt - internal helper function which implements both normMin
##   and normMax to avoid code duplication.
##
normOpt <- function(LSL, USL, maximize=FALSE) {
  if (LSL >= USL) stop("LSL must be less than USL.")
  if (!is.numeric(LSL)) stop("LSL must be numeric.")
  if (!is.numeric(USL)) stop("USL must be numeric.")
  
  ev <- if (maximize) {
    function(y, ...) pnorm((y-a)/b)
  } else {
    function(y, ...) 1 - pnorm((y-a)/b)
  }
  a <- (LSL + USL) / 2
  b <- (USL - a) / 2
  attr(ev, "y.range") <- c(LSL, USL)
  attr(ev, "desire.type") <- ifelse(maximize, "Maximizing Normal", "Minimizing Normal")
  class(ev) <- c(ifelse(maximize, "normMax", "normMin"), "desire.function")
  rm(LSL, USL)
  return(ev)
}

normMax <- function(LSL, USL) 
  normOpt(LSL, USL, maximize=TRUE)

normMin <- function(LSL, USL) 
  normOpt(LSL, USL, maximize=FALSE)

normTarget <- function(LSL, T, USL) {
  if (LSL >= T) stop("LSL must be smaller than T.")
  if (T >= USL) stop("T must be smaller than USL.")
  if (!is.numeric(LSL)) stop("LSL must be numeric.")
  if (!is.numeric(T)) stop("T must be numeric.")
  if (!is.numeric(USL)) stop("USL must be numeric.")
  
  ev <- function(y, ...)
    sqrt(dl(y) * dr(y))
  
  dl <- normMax(LSL, T)
  dr <- normMin(T, USL)
  attr(ev, "y.range") <- c(LSL, USL)
  attr(ev, "desire.type") <- "Target Normal"
  class(ev) <- c("normTarget","desire.function")
  return(ev)
}

ddesire.normMin <- function(x, f, mean=0, sd=1) {
  fun <- function(u, a, b) {
    if (u <= 0 || u >= 0.977249868051820791) {
      return(0)
    } else {
      c1 <- dnorm(qnorm(1-u), mean=(mean - a)/b, sd=sd/b)
      c2 <- dnorm(qnorm(1-u))
      return(c1/c2)
    }    
  }
  env <- environment(f)
  sapply(x, fun, a=env$a, b=env$b)
}

ddesire.normMax <- function(x, f, mean=0, sd=1) {
  fun <- function(u, a, b) {
    if (u <= 0 || u >= 0.977249868051820791) {
      return(0)
    } else {
      c1 <- dnorm(qnorm(u), mean=(mean - a)/b, sd=sd/b)
      c2 <- dnorm(qnorm(u))
      return(c1/c2)
    }    
  }
  env <- environment(f)
  sapply(x, fun, a=env$a, b=env$b)
}

##
## ddesire.normTarget
##
## OME FIXME: Numerically unstable!
##
ddesire.normTarget <- function(x, f, mean=0, sd=1) {
  fun <- function(u, lower, upper) {
    if (u <= 0 || u >= 0.9722) 
      return(0)

    if (u < f(lower) && f(lower) < f(upper)) {
      y <- lower
    } else if (u < f(upper) && f(upper) < f(lower)) {
      y <- upper
    } else {
      y <- uniroot(function(z) f(z) - u, c(lower, upper))$root
    }

    phil <- dnorm(y, al, bl)
    phir <- dnorm(y, ar, br)
    Phil <- pnorm(y, al, bl)
    Phir <- pnorm(y, ar, br)
    
    c1 <- phil/sqrt(Phil)*sqrt(1-Phir)
    c2 <- phir/sqrt(1-Phir)*sqrt(Phil)
    
    res <- dnorm(y, mean, sd) * 2/(c1 - c2)
    return(res)
  }
  env <- environment(f)
  LSL <- env$LSL
  T <- env$T
  USL <- env$USL
  al <- environment(env$dl)$a
  bl <- environment(env$dl)$b
  ar <- environment(env$dr)$a
  br <- environment(env$dr)$b
  l <- sapply(x, fun, lower=9*LSL, upper=T)
  r <- sapply(x, fun, lower=T, upper=9*USL)
  return(l-r)
}

pdesire.normTarget <- function(q, f, mean, sd) {
  fun <- function(u) {
    yl <- uniroot(function(z) f(z) - u, c(-Inf, T))
    yr <- uniroot(function(z) f(z) - u, c(T, Inf))
    if (u <= 0 || u >= 0.9722) {
      return(0)
    } else {
      l <- pnorm(yl, mean, sd)
      r <- 1 - pnorm(yr, mean, sd)
      return(l + r)
    }
  }
  env <- environment(f)
  LSL <- env$LSL
  USL <- env$USL
  T <- env$T
  al <- environment(env$dl)$a
  bl <- environment(env$dl)$b
  ar <- environment(env$dr)$a
  br <- environment(env$dr)$b
  p <- sapply(q, fun)
  return(p)
}

edesire.normMin <- function(f, mean, sd) {
  fun <- function(x)
    x * ddesire(x, f, mean, sd)
  integrate(fun, 0, 0.9722, abs.tol=0)$value
}

edesire.normMax <- edesire.normMin
