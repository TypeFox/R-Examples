#' @title
#' Exponent measure of the HKEVP
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' @description 
#' Computation of the exponent measure \eqn{V(z_1,...,z_n)} of the HKEVP of \cite{Reich and Shaby (2012)}, with given model parameters.
#' 
#' 
#' @param z
#' A numerical vector.
#' Vector \eqn{(z_1,...,z_n)} where the exponent measure is computed.
#' 
#' @param sites
#' A matrix of real values.
#' Coordinates of the sites where the exponent measure is evaluated. Each row corresponds to a position and each column to a coordinate.
#' 
#' @param knots
#' A numerical matrix of real values.
#' Coordinates of the knots: each row corresponds to a knot position and each column to a spatial coordinate.
#' 
#' @param alpha
#' A single numerical value in (0,1].
#' Dependence parameter \eqn{\alpha} in the HKEVP. Low value (resp. value close to 1) corresponds to the limit case of total spatial dependence (resp. independence).
#' 
#' @param tau
#' A positive single numerical value.
#' Bandwidth parameter \eqn{\tau} of the kernel functions in the HKEVP.
#' 
#' 
#' 
#' 
#' 
#' @details
#' The exponent measure describes the spatial dependence structure of a max-stable process, independently from the values of the marginal parameters. If \eqn{Z(\cdot)} is a simple max-stable process, i.e. with unit GEV(1,1,1) margins, recorded at the set of sites \eqn{(s_1, \ldots, s_n)}, its joint cumulative probability density function is given by:
#' \deqn{P\{ Z(s_1)<z_1, \ldots, Z(s_n)<z_n \} = \exp(-V(z_1, \ldots, z_n)) ~,}
#' where \eqn{V} is the so-called exponent measure.

#' For the HKEVP, the exponent measure is explicit for any number \eqn{n} of sites :
#' \deqn{V(z_1, \ldots, z_n) = \sum_{\ell=1}^L \left[ \sum_{i=1}^n \left(\frac{\omega_\ell(s_i)}{z_i}\right)^{1/\alpha}\right]^{\alpha} ~.}
#' 
#' 
#'
#' @return
#' A numerical single value.
#' 
#' 
#' @seealso \link{posterior.expmeasure}
#' 
#' 
#' @export
#' 
#' 
#' @references 
#' Reich, B. J., & Shaby, B. A. (2012). A hierarchical max-stable spatial model for extreme precipitation. The annals of applied statistics, 6(4), 1430. <DOI:10.1214/12-AOAS591>
#' 
#'
#' 
#' 
#' 
#' @examples
#' sites <- as.matrix(expand.grid(1:3, 1:3))
#' knots <- sites
#' alpha <- .4
#' tau <- 1
#' z <- rep(1,9)
#' hkevp.expmeasure(z, sites, knots, alpha, tau)
#' 
#' 
#' 
#' 
hkevp.expmeasure <- function(z, sites, knots, alpha, tau) {
  # Catch Errors and default values
  if (alpha <= 0 | alpha > 1) stop("alpha must be between 0 and 1!")
  if (tau <= 0) stop("tau must be positive!")
  if (length(z) != nrow(sites)) stop("z and sites must have same length!")
  
  # Computing the kernel matrix omega
  nsites <- nrow(sites)
  dsk <- as.matrix(dist(rbind(sites, knots)))[1:nsites, -(1:nsites)]
  omega <- exp(-dsk^2/(2*tau^2))
  omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
  
  # Computing the exponent measure
  result <- sum(colSums(sweep(x = omega, MARGIN = 1, STATS = z, FUN = '/') ^ (1/alpha) ) ^ alpha)
  return(result)
}
