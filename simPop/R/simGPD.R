# Taken from the now archived R package POT
#Descripton of the last available version:
#Package: POT
#Version: 1.1-3
#Date: 2012-10-30
#Title: Generalized Pareto Distribution and Peaks Over Threshold
#Author: Mathieu Ribatet <mathieu.ribatet@math.univ-montp2.fr>
#    Maintainer: Mathieu Ribatet <mathieu.ribatet@math.univ-montp2.fr>
#    Depends: R (>= 1.8.0)
#Description: Some functions useful to perform a Peak Over Threshold
#analysis in univariate and bivariate cases. A user's guide is
#    available.
#    License: GPL (>= 2)
#    URL: http://r-forge.r-project.org/projects/pot/
#    Repository: CRAN
#    Repository/R-Forge/Project: pot
#    Repository/R-Forge/Revision: 492
#    Repository/R-Forge/DateTimeStamp: 2012-10-30 14:21:03
#    Date/Publication: 2012-11-06 09:49:26
#    Packaged: 2012-10-30 15:22:42 UTC; rforge

## This file contains several functions to simulate pseudo-random
## numbers GPD distributed, evaluate the Cumulative Distribution
## Function, evaluate the Quantile Function and the density.

rgpd <- function(n, loc = 0, scale = 1, shape = 0){
  
    if (min(scale) < 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if (shape == 0) 
        return(loc + scale * rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1)/shape)

  }

qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  lambda = 0){
  
    if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
        1) 
      stop("`p' must contain probabilities in (0,1)")
    if (min(scale) < 0) 
      stop("invalid scale")
    if (length(shape) != 1) 
      stop("invalid shape")
    if ((lambda < 0) || (lambda >= 1) || length(lambda) != 1)
      stop("invalid lambda")
    if (any(p < lambda))
      stop("``p'' must satisfy ``p >= lambda''")
    if (lower.tail) 
        p <- 1 - p
    p <- p / (1 - lambda)
    if (shape == 0) 
      return(loc - scale * log(p))
    else return(loc + scale * (p^(-shape) - 1)/shape)

  }

dgpd <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE){
  
    if (min(scale) <= 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    d <- (x - loc)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if (shape == 0) {
        d[index] <- log(1/scale[index]) - d[index]
        d[!index] <- -Inf
    }
    else {
        d[index] <- log(1/scale[index]) - (1/shape + 1) * log(1 + 
            shape * d[index])
        d[!index] <- -Inf
    }
    if (!log) 
        d <- exp(d)

    return(d)

  }

pgpd <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  lambda = 0){
  
    if (min(scale) <= 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if ((lambda < 0) || (lambda >= 1) || length(lambda) != 1)
      stop("invalid lambda")
    q <- pmax(q - loc, 0)/scale
    if (shape == 0) 
        p <- 1 - (1 - lambda) * exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - (1 - lambda) * p^(-1/shape)
    }
    if (!lower.tail) 
        p <- 1 - p

    return(p)

  }





