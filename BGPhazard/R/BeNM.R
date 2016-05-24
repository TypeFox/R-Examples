BeNM <-
function(times, delta) {
  tao <- BeTao(times)
  K <- max(times)
  t.unc <- sort(times[delta == 1])
  n <- rep(0,  K)
  for(j in 1:length(t.unc)) {
    k <- 1
    while (tao[k] != t.unc[j]) {
      k <- k + 1
    }
    n[k] <- n[k] + 1
  }
  m <- rep(0,  K)
  for(k in 1:K) {
    m[k] <- length(times[times > tao[k]])
  }
  m[K] <- 0
  out <- list(n = n, m = m, tao = tao, K = K, t.unc = t.unc)
  out
}
