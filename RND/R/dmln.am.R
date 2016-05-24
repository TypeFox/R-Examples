dmln.am <-
function(x, u.1, u.2, u.3, sigma.1, sigma.2, sigma.3, p.1, p.2)
{
  p.3    = 1 - p.1 - p.2
  out   = p.1 * dlnorm(x, meanlog = u.1, sdlog = sigma.1, log = FALSE) + 
          p.2 * dlnorm(x, meanlog = u.2, sdlog = sigma.2, log = FALSE) + 
          p.3 * dlnorm(x, meanlog = u.3, sdlog = sigma.3, log = FALSE) 
  out

}
