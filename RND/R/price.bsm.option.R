price.bsm.option <-
function(s0, k, r, te, sigma, y)
{
  d1 = (log(s0/k) + (r - y + (sigma^2)/2) * te) / (sigma * sqrt(te))
  d2 = d1 - sigma * sqrt(te)
  cp = s0 * exp(-y*te) * pnorm(d1) - k * exp(-r*te) * pnorm(d2)
  pp = k * exp(-r*te) * pnorm(-d2) - s0 * exp(-y*te) * pnorm(-d1)
  out = list(d1 = d1, d2 = d2, call = cp, put = pp)
  out
}
