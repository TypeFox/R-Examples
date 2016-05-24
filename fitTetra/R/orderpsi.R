orderpsi <-
function(psi, o) {
  psi$mu <- psi$mu[o]
  psi$sigma <- psi$sigma[o]
  psi$p <- psi$p[o]
  psi
}
