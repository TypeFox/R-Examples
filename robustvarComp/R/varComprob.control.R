#############################################################
# 
#	varComprob.control function
#	Author: Claudio Agostinelli and Victor J. Yohai
#	E-mail: claudio@unive.it
#	Date: July, 01, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#                      and Victor J. Yohai
#
#############################################################

varComprob.control <- function(init=NULL, lower=0, upper=Inf, epsilon = 0.001,
    tuning.chi = NULL, bb = 0.5, tuning.psi = NULL,
    arp.chi = 0.1, arp.psi = NULL, max.it = 100,
    rel.tol.beta = 1e-6, rel.tol.gamma = 1e-5, rel.tol.scale = 1e-5,
    trace.lev = 0,
    method = c('compositeTau', "compositeS", "compositeMM", "Tau", "S", "MM"), 
    psi = c('optimal', 'bisquare', 'rocke'),
    beta.univ = FALSE, gamma.univ = FALSE,
    fixed.init=c("lmrob.S", "lmRob"),
    cov.init=c('2SGS', 'covOGK'), cov=TRUE, ...) {

  method <- match.arg(method)
  psi <- match.arg(psi)
  fixed.init <- match.arg(fixed.init)
  cov.init <- match.arg(cov.init)

  if (is.null(init))
    init <- list()
##  if (any(method==c("compositeTau", "compositeMM", "compositeS") & psi=="rocke"))
##    stop("Rocke rho function is not yet available for composite methods")
  if (missing(tuning.chi) || is.null(tuning.chi))
    tuning.chi <- switch(psi,
      'bisquare' = 2.66,
      'rocke' = NULL,
      'optimal' = 1)
  if (missing(tuning.psi) || is.null(tuning.psi))
    tuning.psi <- switch(psi,
      'bisquare' = 3.51,
      'rocke' = NULL,
      'optimal' = sqrt(2.7))
  
  c(list(init = init, lower = lower, upper = upper, epsilon=epsilon,
    psi = psi, tuning.chi = tuning.chi, bb = bb, tuning.psi = tuning.psi,
    arp.chi = arp.chi, arp.psi = arp.psi, max.it = max.it,
    rel.tol.beta = rel.tol.beta, rel.tol.gamma = rel.tol.gamma,
    rel.tol.scale = rel.tol.scale,
    trace.lev = trace.lev, method = method, cov = cov,
    fixed.init=fixed.init, cov.init=cov.init,
    beta.univ = beta.univ, gamma.univ = gamma.univ), list(...))
}
