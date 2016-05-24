compute.implied.volatility <-
function(r, te, s0, k, y, call.price, lower, upper)
{
  n = length(k)
  sigma = numeric(n)
  for (i in 1:n)
  {
    f = function(x) 
    {
         d1 = (log(s0/k[i]) + (r - y + (x^2)/2) * te) / (x * sqrt(te))
         d2 = d1 - x * sqrt(te)
        out = s0 * exp(-y*te) * pnorm(d1) - k[i] * exp(-r*te) * pnorm(d2) - call.price[i]
        out
    }

  sigma[i] = uniroot(f, lower = lower, upper = upper, maxiter = 5000)$root
  }
  sigma 
}
