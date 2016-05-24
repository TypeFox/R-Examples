##============================================================================
## These functions are adapataions of [r,d,p,q]gpd fucntions from the
## evd package by Stephenson et al.
##
## The main difference with evd::functions is that these functions
## return NaNs when the scale parameter is negative, and thus can be
## used in non-contrained opiimisation. 
## 
## Note also that the quantile function is alowed to have p equal to
## 0 or 1, returning then the end-points of the distribution.
## ============================================================================

rGPD <- function(n, loc = 0.0, scale = 1.0, shape = 0.0) {

    if (length(scale) != 1) stop("'scale' must have length 1")
    if (length(shape) != 1) stop("'shape' must have length 1")
    if (scale <= 0) return(rep(NaN, n))
    
    if (shape == 0) return(loc + scale*rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1) / shape)
}

dGPD <- function(x, loc = 0.0, scale = 1.0, shape = 0.0, log = FALSE) {
    
    if (length(scale) != 1) stop("'scale' must have length 1")
    if (length(shape) != 1) stop("invalid shape")
    if (scale <= 0) return(rep(NaN, length(x)))
    
    d <- (x - loc) / scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)

    if (shape == 0) {
      d[index] <- log(1/scale[index]) - d[index]
      d[!index] <- -Inf
    } else {
        d[index] <- log(1/scale[index]) - (1/shape + 1) *
          log(1 + shape * d[index])
        d[!index] <- -Inf
    }
    if (!log) d <- exp(d)
    d
}

pGPD <- function(q, loc = 0.0, scale = 1.0, shape = 0.0, lower.tail = TRUE){
    
    if (length(scale) != 1) stop("'scale' must have length 1")
    if (length(shape) != 1) stop("invalid shape")
    if (scale <= 0) return(rep(NaN, length(q)))
    
    q <- pmax(q - loc, 0) / scale
    if (shape == 0) p <- 1 - exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - p^(-1 / shape)
    }
    if(!lower.tail) p <- 1 - p
    p
}

qGPD <-function(p, loc = 0.0, scale = 1.0, shape = 0.0, lower.tail = TRUE) {

    if (length(scale) != 1L) stop("'scale' must have length 1")
    if (length(shape) != 1L) stop("invalid shape")
    if (scale <= 0) return(rep(NaN, length(p)))
     
    if (min(p, na.rm = TRUE) < 0 || max(p, na.rm = TRUE) > 1)
        stop("`p' must contain probabilities in (0, 1)")

    if (lower.tail) p <- 1 - p
    if (shape == 0) {
        return(loc - scale * log(p))
    } else {
        return(loc + scale * (p^(-shape) - 1) / shape)
    }
}

HGPD <- function(x, loc = 0.0, scale = 1.0, shape = 0.0) {
    if (length(scale) != 1L) stop("'scale' must have length 1")
    if (length(shape) != 1L) stop("invalid shape")
    if (scale <= 0) return(rep(NaN, length(x)))
    
    -log(pGPD(x, loc = loc, scale = scale, shape = shape,
              lower.tail = FALSE))
}

hGPD <- function(x, loc = 0.0, scale = 1.0, shape = 0.0) {
    if (length(scale) != 1L) stop("'scale' must have length 1")
    if (length(shape) != 1L) stop("invalid shape")
    if (scale <= 0) return(rep(NaN, length(x)))

    h <- rep(0, length(x))
    if (shape == 0) {
        index <- (x >= loc)
        if (any(index)) h[index] <- 1 / scale
    } else if (shape < 0) {
        index <- (x >= loc) & (x <= loc - scale / shape) 
        if (any(index)) h[index] <- 1 / (scale + shape * (x[index] - loc))
    } else {
        index <- (x >= loc) 
        if (any(index)) h[index] <- 1 / (scale + shape * (x[index] - loc))
    }
    h
}
