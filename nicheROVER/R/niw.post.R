niw.post <-
function(nsamples, X, lambda, kappa, Psi, nu) {
  par <- niw.coeffs(X, lambda, kappa, Psi, nu)
  rniw(nsamples, par$lambda, par$kappa, par$Psi, par$nu)
}
