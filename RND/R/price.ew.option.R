price.ew.option <-
function(r, te, s0, k, sigma, y, skew, kurt)
{
  discount.factor  = exp(-r * te)

  v                =  sqrt(exp(sigma^2 * te) - 1)
  m                =  log(s0) + (r - y - (sigma^2)/2) * te
  skew.lognorm     =  3 * v + v^3
  kurt.lognorm     =  16 * v^2 + 15 * v^4 + 6 * v^6 + v^8
  cumul.lognorm    =  (s0 * exp((r-y)*te) * v)^2
  
  density.lognorm  =  dlnorm(x = k, meanlog = log(s0) + (r - y - (sigma^2)/2)*te, sdlog = sigma * sqrt(te), log = FALSE)
  frst.der.lognorm =  -1 * ( 1 + (log(k) - m)/(te * sigma^2) ) * density.lognorm/k
  scnd.der.lognorm =  -1 * ( 2 + (log(k) - m)/(te * sigma^2) ) * frst.der.lognorm / k  - density.lognorm/(k^2 * sigma^2)
  call.lognorm     =  price.bsm.option(r = r, te = te, s0 = s0, k = k, sigma = sigma, y = y)$call

  ew.call =  call.lognorm - discount.factor * (skew - skew.lognorm) * cumul.lognorm^(1.5) * frst.der.lognorm/6 +
                                      discount.factor * (kurt - kurt.lognorm) * cumul.lognorm^2 * scnd.der.lognorm/24
  ew.put  =  ew.call + k * exp(-r * te) - s0 * exp(-y * te)

  out = list(call = ew.call, put = ew.put)
  out
   
}
