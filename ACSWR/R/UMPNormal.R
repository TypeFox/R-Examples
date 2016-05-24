UMPNormal <-
function(mu0, sigma, n,alpha)  {
  mu0-qnorm(alpha)*sigma/sqrt(n)
}
