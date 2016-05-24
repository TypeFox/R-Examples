dmln <-
function(x, alpha.1, meanlog.1, meanlog.2, sdlog.1, sdlog.2)
{
  alpha.2 = 1 - alpha.1
  out   = alpha.1 * dlnorm(x, meanlog = meanlog.1, sdlog = sdlog.1, log = FALSE) + 
          alpha.2 * dlnorm(x, meanlog = meanlog.2, sdlog = sdlog.2, log = FALSE)
  out

}
