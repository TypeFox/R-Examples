#' @title
#' Simulation of the HKEVP
#' 
#' @author
#' Quentin Sebille
#'
#'
#'
#' @description
#' Simulation procedure of the HKEVP with given sites and knots positions and marginal and spatial dependence parameters.
#'
#'
#'
#'
#'
#' @param nrep
#' A positive integer.
#' Number of realisations of the block maxima process.
#' 
#' @param sites
#' A numerical matrix of real values.
#' Coordinates for the sites positions: each row corresponds to a position and each column to a spatial coordinate.
#' 
#' @param knots
#' A numerical matrix of real values.
#' Coordinates of the knots: each row corresponds to a knot position and each column to a spatial coordinate. If missing, \code{knots} is equal to \code{sites}.
#' 
#' @param mu
#' A numerical value or a vector of real values.
#' GEV location parameter. If \code{mu} is a vector, its length must coincide with the number of sites.
#' 
#' @param sigma
#' A numerical value or a vector of positive values.
#' GEV scale parameter. If \code{sigma} is a vector, its length must coincide with the number of sites.
#' 
#' @param xi
#' A numerical value or a vector of real values.
#' GEV shape parameter. If \code{xi} is a vector, its length must coincide with the number of sites.
#' 
#' @param alpha
#' A numerical value in (0,1].
#' Dependence parameter \eqn{\alpha} in the HKEVP. Low value (resp. value close to 1) corresponds to the limit case of complete spatial dependence (resp. independence).
#' 
#' @param tau
#' A positive numerical value.
#' Bandwidth parameter \eqn{\tau} of the kernel functions in the HKEVP.
#' 
#' 
#' 
#' 
#' 
#' @details
#' Simulating one realisation of the block maxima process \eqn{Y(\cdot)} from the HKEVP involves three steps:
#' \enumerate{
#' \item The \emph{nugget process} \eqn{U(\cdot)} is generated independently at each position, by simulating a random variable with \eqn{GEV(1,\alpha,\alpha)} distribution.
#' \item The \emph{residual dependence process} \eqn{\theta(\cdot)} is computed by using the kernel functions centered at the set of knots, the bandwidth parameter \eqn{\tau} and the simulations of the positive stable \eqn{PS(\alpha)} random effect \eqn{A}.
#' \item The process \eqn{Z = U\theta} is computed and its margins are transformed to the general GEV distribution with \eqn{\mu(s),\sigma(s)} and \eqn{\xi(s)} parameters.
#' }
#' 
#' 
#' 
#' 
#'
#' @return
#' A numerical matrix of real values.
#' Each column corresponds to a position and each row to a realisation of the process.
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
#' 
#' 
#' 
#' 
hkevp.rand <- function(nrep, sites, knots, mu, sigma, xi, alpha, tau) {
  # Catching errors
  if (!is.matrix(sites)) stop("sites must be a matrix!")
  if (missing(knots)) knots <- sites
  if (!is.matrix(knots)) stop("knots must be a matrix!")
  if (length(alpha) != 1) stop("alpha must be of length one!")
  if (alpha <= 0 | alpha > 1) stop("alpha must be within (0,1] !")
  if (length(tau) != 1) stop("tau must be of length one!")
  if (tau <= 0) stop("tau must be positive!")
  if (any(sigma <= 0)) stop("sigma must be positive!")
  if (ncol(sites) != ncol(knots)) stop("sites and knots must have same dimensions!")
  
  # Number of sites and knots
  nsites <- nrow(sites)
  nknots <- nrow(knots)
  dsk <- as.matrix(dist(rbind(sites, knots)))[1:nsites, -(1:nsites)]
  
  # Positive Stable random variable generation
  unifPS <- matrix(runif(nrep*nknots, 0, pi), nrep, nknots)
  expoPS <- matrix(rexp(nrep*nknots, 1), nrep, nknots)
  A1 <- (sin((1 - alpha)*unifPS) / expoPS) ^ ((1 - alpha)/alpha)
  A2 <- sin(alpha*unifPS) / (sin(unifPS)) ^ (1/alpha)
  A <- A1*A2
  
  # Theta and U processes
  omega <- exp(-dsk^2/(2*tau^2))
  omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
  theta <- (A %*% t(omega^(1/alpha)))^alpha
  U <- matrix(rexp(nrep*nsites) ^ (-alpha), nrep, nsites) # ~ GEV(1,alpha,alpha)
  Z <- theta*U
  
  # Transformation into GEV margins
  gevloc <- matrix(mu, nsites, 1)
  gevscale <- matrix(sigma, nsites, 1)
  gevshape <- matrix(xi, nsites, 1)
  result <- Z
  for (i in 1:nsites) {
    if (gevshape[i] == 0)
      result[,i] <- gevloc[i] + gevscale[i] * log(Z[,i])
    else
      result[,i] <- gevloc[i] + (gevscale[i]/gevshape[i]) * (Z[,i] ^ (gevshape[i]) - 1)
  }
  
  return(result)
}

