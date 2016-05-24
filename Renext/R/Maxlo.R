##======================================================
## The name of this distibution is a mnemonic reference
## to Lomax
##======================================================

rmaxlo <- function(n, scale = 1.0, shape = 4.0) {
  
  ## if (scale <= 0) stop("'scale' must be > 0")
  ## if (shape <= 0) stop("'shape' must be > 0")
  if ((scale <= 0) || (shape <= 0)) return(rep(NaN, n))
  X <- rexp(n, rate = shape)
  X <- scale*(1-exp(-X))
}

dmaxlo <- function(x, scale = 1.0, shape = 4.0, log = FALSE) {
  
  ## if (scale <= 0) stop("'scale' must be > 0")
  ## if (shape <= 0) stop("'shape' must be > 0")
  if ((scale <= 0) || (shape <= 0)) return(rep(NaN, length(x)))
  f <- rep(0, length(x))
  ind <- (x > 0) & (x < scale)
  f[ind] <- log(shape/scale) + (shape - 1.0)*log(1 - x[ind]/scale)
  if (!log) f[ind] <- exp(f[ind])
  f
  
}

pmaxlo <- function(q, scale = 1.0, shape = 4.0, lower.tail = TRUE) {

  ## if (scale <= 0) stop("'scale' must be > 0")
  ## if (shape <= 0) stop("'shape' must be > 0")
  if ((scale <= 0) || (shape <= 0)) return(rep(NaN, length(q)))
  F <- rep(0, length(q))
  ind <- (q > 0) 
  F[ind] <- (1 - q[ind] / scale)^shape
  if (lower.tail) F[ind] <- 1 - F[ind]
  else F[!ind] <- 1
  F
  
}

qmaxlo <- function(p, scale = 1.0, shape = 4.0) {

  ## if (scale <= 0) stop("'scale' must be > 0")
  ## if (shape <= 0) stop("'shape' must be > 0")
  if ((scale <= 0) || (shape <= 0)) return(rep(NaN, length(p)))
  if ( any(p < 0) || any(p > 1) ) stop("'p' must be >=0 and <=1")
  scale* (1 - (1-p)^(1/shape))

}
