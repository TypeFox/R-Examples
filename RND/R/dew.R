dew <-
function(x, r, y, te, s0, sigma, skew, kurt)
{
  v                =  sqrt(exp(sigma^2 * te) - 1)
  m                =  log(s0) + (r - y - (sigma^2)/2) * te
  skew.lognorm     =  3 * v + v^3
  kurt.lognorm     =  16 * v^2 + 15 * v^4 + 6 * v^6 + v^8
  cumul.lognorm    =  (s0 * exp((r-y) * te) * v)^2

  density.lognorm  =  dlnorm(x = x, meanlog = log(s0) + (r - y - (sigma^2)/2)*te, sdlog = sigma * sqrt(te), log = FALSE)
  frst.der.lognorm =  -1 * ( 1 + (log(x) - m)/(te * sigma^2) ) * density.lognorm/x
  scnd.der.lognorm =  -1 * ( 2 + (log(x) - m)/(te * sigma^2) ) * frst.der.lognorm / x  - density.lognorm/(x^2 * sigma^2)
  thrd.der.lognorm =  -1 * ( 3 + (log(x) - m)/(te * sigma^2) ) * scnd.der.lognorm / x  - 2 * frst.der.lognorm/(x^2 * sigma^2) + density.lognorm/(x^3 * sigma^2)
  frth.der.lognorm =  -1 * ( 4 + (log(x) - m)/(te * sigma^2) ) * thrd.der.lognorm / x  - 3 * scnd.der.lognorm/(x^2 * sigma^2) + 
                            3 * frst.der.lognorm/(x^3 * sigma^2) - 2 * density.lognorm/(x^4 * sigma^2)
  
  out = density.lognorm - (skew - skew.lognorm) * ((cumul.lognorm)^(1.5))*thrd.der.lognorm/6 + (kurt - kurt.lognorm)*((cumul.lognorm)^2)*frth.der.lognorm/24

}
