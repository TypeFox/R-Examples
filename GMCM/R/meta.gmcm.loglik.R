# A wrapper for dgmcm.loglik in the meta analysis case
# where the loglik have arguments in the transformed parameter space.
meta.gmcm.loglik <- function (tpar, u, rescale = TRUE, positive.rho = TRUE) {
  # par is a vector of transformed (if rescale == TRUE) parameters
  if (rescale) {
    par <- tt(tpar, ncol(u), positive.rho = positive.rho)
  } else {
    par <- tpar
  }
  return(dgmcm.loglik(theta = meta2full(par, ncol(u)), u))
}
