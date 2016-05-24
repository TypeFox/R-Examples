### Generate random phi values from lognormal distribution and matched by the
### rank of SCUO indices.

scuo.random <- function(SCUO, phi.Obs = NULL,
    meanlog = .CF.PARAM$phi.meanlog, sdlog = .CF.PARAM$phi.sdlog){
  if(!is.null(phi.Obs)){
    meanlog <- mean(log(phi.Obs))
    sdlog <- sqrt(var(log(phi.Obs)))
  }

  n <- length(SCUO)
  x <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  ret <- sort(x)[rank(SCUO)]

  ret
} # End of scuo.random().
