price.shimko.option <-
function(r, te, s0, k, y , a0, a1, a2)
{
  
  sigma.k = a0 + a1 * k + a2 * k^2
  mu      = log(s0) + (r - y - (sigma.k^2)/2) * te
  integrand = function(x) { (1/x)*(x - k) * (1/(sqrt(2*pi*te)*sigma.k))*exp( -1*((log(x) - mu)^2)/(2 * (sigma.k^2) * te ) ) }
  call = exp(-r * te) * integrate(integrand, k, Inf, subdivisions=500)$value
  put  = call + exp(-r * te) * k - exp(-y*te)* s0
  out  = list(call = call, put = put)
  out

}
