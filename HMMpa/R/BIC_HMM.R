BIC_HMM <-
function(size, m, k, logL)  
{ 
  foo <- try(-2 * logL + log(size) * (m^2 + k * m - 1), silent = FALSE)
  if (inherits(foo, "try-error")) 
  {
    foo <- NA
  }
  return(foo)
}
