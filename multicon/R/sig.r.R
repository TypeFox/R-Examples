sig.r <-
function(r, n, tail) {
  crit.p <- c(.05, .025, .005, .0005)
  if (tail == 1) crit.p <- crit.p * 2
  crit.t <- qt(crit.p, n-2, lower.tail=FALSE) # assumes that n is a single number
  crit.r <- sqrt( crit.t^2 / (crit.t^2 + (n-2)))
  cut(r, c(-1, crit.r, 1.01), labels=c("   ", "+  ", "*  ", "** ", "***"))
}
