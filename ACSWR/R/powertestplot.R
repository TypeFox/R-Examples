powertestplot <-
function(mu0,sigma,n,alpha)  {
  mu0seq <- seq(mu0-3*sigma, mu0+3*sigma,(6*sigma/100))
  betamu <- pnorm(sqrt(n)*(mu0-mu0seq)/sigma-qnorm(1-alpha))
  plot(mu0seq,betamu,"l",xlab=expression(mu),ylab="Power of UMP Test",main=expression(paste("H:",mu >= mu[0]," vs K:",mu<mu[0])))
  abline(h=alpha)
  abline(v=mu0)
}
