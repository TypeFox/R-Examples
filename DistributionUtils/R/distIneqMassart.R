distIneqMassart <- function(densFn = "norm", n = 10000,
                            probBound = 0.001, ...){
  ## Purpose: Sample from distribution and carry out Massart test
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott and Christine Yang Dong, Date:  8 Feb 2010, 22:07

  ## dpqr random test based on the test in base R
  ## RNG tests using DKW inequality for rate of convergence
  ##
  ## P(sup | F_n - F | > t) < 2 exp(-2nt^2)
  ##
  ## The 2 in front of exp() was derived by Massart. It is the best possible
  ## constant valid uniformly in t,n,F. For large n*t^2 this agrees with the
  ## large-sample approximation to the Kolmogorov-Smirnov statistic.
  ##

  dfun <- match.fun(paste("d", densFn, sep = ""))
  pfun <- match.fun(paste("p", densFn, sep = ""))
  rfun <- match.fun(paste("r", densFn, sep = ""))

  ## sample from distribution
  x <- rfun(n = n, ...)

  ## calculate sup of empirical from true distribution function
  tx <- table(x)
  xi <- as.numeric(names(tx))
  f <- pfun(xi, ...)
  fhat <- cumsum(tx)/n
  sup <- max(abs(fhat - f))
  supError <- signif(sup, 2)

  ## probability of that sup value
  pVal <- min(1, round(2*exp(-2*n*sup*sup), 4))

  ## value of t for that probability bound
  tVal <- sqrt(log(probBound/2)/(-2*n))

  ## check if inequality satisfied
  if (sup < sqrt(log(probBound/2)/(-2*n))){
    check <- TRUE
  } else {
    check <- FALSE
  }

  results <- list(sup = supError, probBound = probBound, t = tVal,
                  pVal = pVal, check = check)

  return(results)
}
