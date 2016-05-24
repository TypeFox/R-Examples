CGaM <-
function(times, delta, type.t, K, covar, theta) {
  tao <- GaTao(times, delta, type.t, K)
  if (type.t == 2) {
    K <- ceiling(max(times))
  }
  m <- rep(0, K)
  for(i in 1:length(times)) {
    for(k in 1:K) {
      if (tao[k + 1] < times[i]) {
        m[k] <- m[k] + (tao[k + 1] - tao[k]) * exp(theta %*% covar[i, ])
      }
      if (tao[k] < times[i] && times[i] <= tao[k + 1]) {
        m[k] <- m[k] + (times[i] - tao[k]) * exp(theta %*% covar[i, ])
      }
    }
  }
  return(m)
}
