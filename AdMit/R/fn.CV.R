## Function which computes the coefficient of variation
## See Hoogerheide (2006, p.48)
## __input__
## w : [Nx1 vector] of weights
## __output__
## [double] coefficient of variation
## __20080429__
'fn.CV' <- function(w)
  {
    r <- sd(w)
    if (r==0)
      stop ("'w' is constant in 'fn.CV'")
    as.numeric(r/mean(w))
  }
