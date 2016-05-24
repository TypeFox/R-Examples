AIC_HMM <- 
function(logL, m, k) 
{ 
  foo <- try(-2 * logL + 2 * (m^2 + k * m - 1), silent = FALSE) 
  if (inherits(foo, "try-error")) 
  {
    foo <- NA
  }
  return(foo)
}
