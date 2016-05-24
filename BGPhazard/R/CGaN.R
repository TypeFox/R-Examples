CGaN <-
function(times, delta, type.t, K, covar) {
  tao <- GaTao(times, delta, type.t, K)
  if (type.t == 2) {
    K <- ceiling(max(times))
  }
  t.unc <- sort(times[delta == 1])
  n <- rep(0, K)
  for(i in 1:length(t.unc)) {
    k <- 1
    while((tao[k] < t.unc[i] && t.unc[i] <= tao[k + 1]) == FALSE) {
      k <- k + 1
    }
    n[k] <- n[k] + 1
  }
  return(n)
}
