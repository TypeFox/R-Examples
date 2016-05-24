#' @title
#' Posterior distribution of the HKEVP exponent measure
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' This function returns the posterior distribution of the exponent measure \eqn{V(z_1, ..., z_n)} evaluated at the set of sites \eqn{(s_1, ..., s_n)}, computed from the output of \code{hkevp.fit}
#'
#'
#'
#'
#'
#' @param fit
#' A named list.
#' Output from the \code{hkevp.fit} procedure.
#' 
#' @param sites
#' A matrix of real values.
#' Coordinates of the sites where the exponent measure is evaluated. Each row corresponds to a position and each column to a coordinate.
#' 
#' @param z
#' A numerical vector.
#' Vector \eqn{(z_1,...,z_n)} where the exponent measure is computed.
#' 
#' 
#' 
#' 
#' 
#' @return
#' A numerical vector.
#' Exponent measure evaluated at each MCMC step.
#' 
#' 
#' @seealso \link{hkevp.expmeasure}
#' 
#' 
#'
#' @examples
#' # Simulation of HKEVP:
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' knots <- sites
#' mu <- sites[,1]*10
#' sigma <- 3
#' xi <- .2
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, sites, knots, mu, sigma, xi, alpha, tau)
#' 
#' # HKEVP fit:
#' fit <- hkevp.fit(ysim, sites, knots, niter = 100, nburn = 50, quiet = FALSE)
#' 
#' ## Exponent measure:
#' z <- rep(1, 2)
#' postem <- posterior.expmeasure(fit, sites[1:2,], z)
#' print(median(postem))
#' print(sd(postem))
#' 
#' 
#' 
#' 
posterior.expmeasure <- function(fit, sites, z) {
  # Catching Error
  if (length(z) != nrow(sites)) stop("'z' and 'sites' does not match!")
  
  nstep <- fit$nstep
  result <- rep(NA, nstep)
  for (i in 1:nstep)
    result[i] <- hkevp.expmeasure(z, sites, fit$knots, fit$alpha[i], fit$tau[i])
  
  return(result)
}