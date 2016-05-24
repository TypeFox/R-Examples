
## -----------------------------------------------------------------------------
## Latin Hypercube distribution
## -----------------------------------------------------------------------------

Latinhyper <- function(parRange, num) {

  npar <- nrow(parRange)
  latin <- NULL
  for (i in 1:npar) {
    pr    <- unlist(parRange[i,])
    ## delta of parameter interval
    dpar  <- diff(pr) / num
    ## index to interval (0=1st interval)
    ii    <- sort(runif(num), index.return = TRUE)$ix - 1
    rr    <- runif(num)
    ## random value within interval
    pval  <- pr[1] + ii*dpar + rr*dpar
    latin <- cbind(latin, pval)
  }
  colnames(latin) <- rownames(parRange)
  latin
}

