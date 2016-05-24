alpha.aci <-
function(x, k, n, CI=.95) {
  pcrit <- (1-CI)/2
  zcrit.L <- qnorm(pcrit,lower.tail=T)
  zcrit.U <- qnorm(pcrit,lower.tail=F)
  w <- sqrt(2*k/(n*(k-1)))
  LL <- 1-(1-x)*exp(zcrit.U*w)
  UL <- 1-(1-x)*exp(zcrit.L*w)
  out <- cbind(LL, UL)
  colnames(out) <- c("Lower Limit", "Upper Limit")
  return(out)
}
