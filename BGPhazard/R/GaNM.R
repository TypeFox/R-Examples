GaNM <-
function(times, delta, type.t, K) {
  tao <- GaTao(times, delta, type.t, K)
  if (type.t == 2) {
    K <- ceiling(max(times))
  }
  t.unc <- sort(times[delta == 1])
  n <- rep(0, K)
  for(j in 1:length(t.unc)) {
    i <- 1
    while((tao[i] < t.unc[j] && t.unc[j] <= tao[i + 1]) == FALSE) {
      i <- i + 1
    }
    n[i] <- n[i] + 1
  }
  m <- rep(0, K)
  for(j in 1:length(times)) {
    for(i in 1:K) {
      if (tao[i + 1] < times[j]) {
        m[i] <- m[i] + tao[i + 1] - tao[i]
      }
      if (tao[i] < times[j] && times[j] <= tao[i + 1]) {
        m[i] <- m[i] + times[j] - tao[i]
      }
    }
  }
  out <- list(n = n, m = m, tao = tao, t.unc = t.unc)
  return(out)
}
