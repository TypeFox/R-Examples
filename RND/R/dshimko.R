dshimko <-
function(r, te, s0, k, y, a0, a1, a2)
{
  sigma = a0 + a1*k + a2*k^2

   v  =  sigma * sqrt(te)
  d1  =  (log(s0/k) + (r - y + (sigma^2)/2) * te) / v
  d2  =  d1 - v
  d1x =  -1/(k * v) + (1 - d1/v)*(a1 + 2 * a2 * k)
  d2x =  d1x - (a1 + 2 * a2 * k)
  out =  -1 * dnorm(d2)*(d2x - (a1 + 2 * a2* k)*(1 - d2*d2x) - 2 * a2 * k)
  out

}
