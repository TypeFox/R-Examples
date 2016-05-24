## -----------------------------------------------------------------------------
## Draws a Uniform Random Sample from a set of parameters
## -----------------------------------------------------------------------------

Unif <- function(parRange, num) {
  npar <- nrow(parRange)
  parset   <- NULL
  for (i in 1:npar)
    parset <- cbind(parset, runif(num, parRange[i, 1], parRange[i, 2]))

  colnames(parset) <- rownames(parRange)
  parset
}

