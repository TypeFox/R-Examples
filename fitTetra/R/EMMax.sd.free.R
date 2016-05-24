EMMax.sd.free <-
function(y,z,mu) {
  sigma <- sqrt(apply(rbind(z,mu),2,wmean2,y))
  sigma[sigma < 0.01] <- 0.01
  sigma
}
