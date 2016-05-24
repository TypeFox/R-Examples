edgeworth <- function(x, n, rho3, rho4, mu, sigma2, deg=3,
                      type=c("standardized", "mean", "sum")) {
  MAX.DEGREE <- 3
  if (!deg %in% 1:MAX.DEGREE) {
    stop("Argument 'degree' must be an integer between 1 and ",
         MAX.DEGREE, "!")
  }
  if (!isInteger(n) || n <=0) {
    stop("Argument 'n' must be a positive integer!")
  }
  if(!is.numeric(sigma2) || sigma2 <= 0) {
    stop("Argument 'sigma2' must be a real value > 0!")
  }
  type <- match.arg(type)
  areMissingMoments <- missing(mu) || missing(sigma2)
  if (areMissingMoments && type != "standardized") {
    stop("For approximating the '", type ,
         "' the first two moments are necessary!")
  }
  areMissingMoments <- (missing(rho3) && deg >= 2) || (missing(rho4) && deg >= 3)
  if (areMissingMoments) {
    stop("Necessary higher moments are missing!")
  }

  switch(type,
         "standardized" = {
           scaling <- 1
           z <- x
         },
         "mean" = {
           scaling <- sqrt(n) / sqrt(sigma2)
           z <- (n * x - n * mu) / sqrt(n * sigma2)
         },
         "sum" = {
           scaling <- 1 / sqrt(n * sigma2)
           z <- (x - n * mu) / sqrt(n * sigma2)
         })
  edgeApprox <- switch(deg,
                       dnorm(z),
                       dnorm(z) * (1 + rho3 / (6 * sqrt(n)) * hermite(z, 3)),
                       dnorm(z) * (1 + rho3 / (6 * sqrt(n)) * hermite(z, 3) +
                                   1 / n * (rho4 / 24 * hermite(z, 4) +
                                            rho3 ^ 2 / 72 * hermite(z, 6)))) 
  value <- scaling * edgeApprox
  return(approximation(y=x, approx=value, n=n, type=type, approx.type="Edgeworth"))
}
