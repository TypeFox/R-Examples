transectHolder <- function(..., family="exponential") {
  transectHolder <- list()
  transectHolder$transects <- list(...)
  distname <- tolower(family)
  transectHolder$family <- family
  transectHolder$parameters <-  fitDistances(transectHolder, distname)
  transectHolder$rng <- switch(distname,
                               "beta" = "rbeta",
                               "chi-squared" = "rchisq",
                               "exponential" = "rexp",
                               "f" = "rf",
                               "gamma" = "rgamma",
                               "log-normal" = "rlnorm",
                               "lognormal" = "rlnorm",
                               "negative binomial" = "rnbinom",
                               "poisson" = "rpois",
                               "weibull" = "rweibull",
                               NULL)
  if (is.null(transectHolder$rng))
    stop("Unsupported distribution")
  class(transectHolder) <- "transectHolder"
  return(transectHolder)
}
