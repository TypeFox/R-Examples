qgpd <-
function(p , sigma, xi, u = 0, lower.tail=TRUE, log.p=FALSE ){
  ## this code is trickier than it looks.

  ## we use G ~ (exp(xi E) - 1) / xi

  ## there's no argument checking in qexp, so we do this to match the
  ## (erroneous) existing behaviour

  if ((!log.p) && any((p < 0) || (p > 1))) {
    stop("p must lie between 0 and 1 if log.p=FALSE")
  }

  ## get the quantiles of the standard exponential
  exp.quantiles <- qexp(p, lower.tail=lower.tail, log.p=log.p)

  ## transform to quantiles of standard gpd
  ## this works for negative xi because we get two sign flips
  standard.gpd  <- .exprel(exp.quantiles * xi) * exp.quantiles

  ## and transform to our gpd.
  u + sigma * standard.gpd
}

