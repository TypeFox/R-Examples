qgev <- function(p, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){
  ## this code is trickier than it looks.

  ## we use G ~ (exp(-xi log E) - 1) / xi

  ## there's no argument checking in qexp, so we do this to match the
  ## existing behaviour

  if ((!log.p) && any((p < 0) || (p > 1))) {
    stop("p must lie between 0 and 1 if log.p=FALSE")
  }

  ## get the quantiles of the standard exponential
  ## change the sense of lower tail (the minus sign above)

  neg.exp.quantiles <- -log(qexp(p, lower.tail=!lower.tail,
                                 log.p=log.p))
  standard.gev <- .exprel(neg.exp.quantiles * xi) * neg.exp.quantiles

  ## and now shift and scale
  mu + sigma * standard.gev
}


