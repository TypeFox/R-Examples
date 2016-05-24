## Function which computes the weigths (exponential)
## __input__
## lnk : [Nx1 vector] of log kernel values
## lnd : [Nx1 vector] of log mixture values
## __output__
## [Nx1 vector] of weights
## __20080429__
'fn.computeexpw' <- function(lnk, lnd)
  {
    r <- lnk-lnd
    r <- r-max(r) ## robustify
    r <- exp(r)
    r[is.na(r) | is.nan(r)] <- 0
    as.vector(r)
  }
