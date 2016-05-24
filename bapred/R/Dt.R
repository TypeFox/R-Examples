Dt <-
function(rho) {
  threshold=0.05
  ut = qnorm(1 - threshold/2)
  delta = unlist(lapply(rho,bivprob,lower=-ut)) - (1 - threshold)^2
  dt = delta/(threshold * (1 - threshold))
  return(dt)
}
