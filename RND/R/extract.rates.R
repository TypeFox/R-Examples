extract.rates <-
function(calls,  puts, s0,  k, te )
{
  y  =  puts - calls
  x  =  k

  linear.model  =  lm(y ~ x)
  slope         =  linear.model$coeff[[2]]
  intercept     =  linear.model$coeff[[1]]

  r = -1 * log(slope) / te
  y = -1 * log( -1 * intercept / s0 ) / te

  out = list(risk.free.rate = r, dividend.yield = y)
  out
}
