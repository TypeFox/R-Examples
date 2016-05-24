chi2cdf <-
function(x, n)
  {
    ###  Rick Aster's matlab function is easy in R
    ###  use pchisq
    return(pchisq(x, n) )

  }
