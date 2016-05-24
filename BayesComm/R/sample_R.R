sample_R <-
function (e, priR) {
  S <- t(e) %*% e + priR[2] * diag(ncol(e))
  v <- priR[1]
  df <- v:(v - ncol(S) + 1)
  sig <- solve(rwish(solve(S), df))
  cov2cor(sig)
}