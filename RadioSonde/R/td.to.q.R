td.to.q <- function(temp, td, p)
{
  # calculate specific humidity (g/kg) using temperature and dew-point T
  es <- 6.1121 * exp((temp * 17.67)/(temp + 243.5))
  e <- 6.1121 * exp((td * 17.67)/(td + 243.5))
  q <- (1000. * (0.62197 * e))/(p - (1. - 0.62197) * e)
  return(q)
}
