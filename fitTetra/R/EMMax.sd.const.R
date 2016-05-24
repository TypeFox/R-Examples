EMMax.sd.const <-
function(y,z,mu) {
  sd <- sqrt( weighted.mean((rep(y,length(mu)) - rep(mu, each=length(y)))^2, as.vector(z)))
  sigma <- rep(sd,length(mu))
  sigma[sigma < 0.01] <- 0.01
  sigma
}
