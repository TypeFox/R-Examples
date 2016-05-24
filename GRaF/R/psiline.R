psiline <-
function(s, adiff, a, K, y, d0, mn = 0, wt) {
  a <- a + s * as.vector(adiff)
  f <- K %*% a + mn
  psi(a, f, mn, y, d0, wt)
}
