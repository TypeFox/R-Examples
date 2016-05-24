## Function which computes the weights
## See Hoogerheide (2006,pp.46-47)
## __input__
## theta  : [Nxk matrix] of draws
## KERNEL : [function] which computes the kernel !! must be vectorized !!
## mit    : [list] containing mixture information
## log     : [logical] natural logartithm output (default: TRUE)
## ...    : additional parameters used by 'KERNEL'
## __output__
## [Nx1 vector] of weights
## __20080427__
'fn.w' <- function(theta, KERNEL, mit, log=TRUE, ...)
  {
    lnk <- KERNEL(theta, log=TRUE, ...)
    lnd <- dMit(theta, mit, log=TRUE)
    if (log)
      {
        r <- lnk-lnd
        r[is.na(r) | is.nan(r)] <- -Inf
      }
    else
      {
        r <- fn.computeexpw(lnk, lnd)
      }
    as.vector(r)
  }
