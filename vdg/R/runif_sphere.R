#' Sampling for hyperspheres/hypercubes
#' 
#' Sample uniformly in or on a hyperspheres or hypercubes. 
#' 
#' @param n number of points to sample
#' @param m number of design factors
#' @param max.dist maximum distance from origin (L-infinity norm/supremum distance) for the hypercuboidal 
#' design region (enveloping hypercube)
#' @param max.radius maximum radius of the hyperspherical design region (enveloping hypersphere)
#' @param at logical indicating whether to sample on concentric hyperspheres/hypercubes or not. With this
#' option \code{n} is distributed proportionally across radii / supremum distances so that the density
#' of samples on each concentric hypercube / hypersphere are uniform across the different hyperspheres / hypercubes..
#' @param nr.dist the number of concentric hypercubes to use in case at is \code{TRUE}
#' @param nr.rad number of concentric hyperspheres to sample on in case of \code{at} being \code{TRUE}
#' @author Pieter C. Schoonees
#' @export
#' @examples
#' 
#' set.seed(1234)
#' samp <- runif_sphere(n = 500, at = TRUE)
#' plot(samp, asp = 1)
runif_sphere <- function (n, m = 2, max.radius = sqrt(m), at = FALSE, nr.rad = 21) {
  if(m > 2){
    bvec <- (m - 1:(m - 2))/2
    bvec <- beta(bvec, bvec)
  } else bvec <- 1
  if(at){
    spec.dist <- seq(from = 0, to = max.radius, length = nr.rad)
    nprop <- gamma(m/2 + 1) * 2^((m - 1) * (m - 2)/2 + 1) *  spec.dist^(m - 1) * 
      bvec / (m * pi^(m / 2) * max.radius^m)
    ns <- ceiling(n * nprop / sum(nprop))
    ns[1] <- 1
    rvec <- rep(spec.dist, ns)
    n <- length(rvec)
  } 
  uvec <- runif(n)
  tmp <- matrix(rnorm(n*m), nrow = n, ncol = m)
  if(!at) rvec <- max.radius * (m * pi^(m/2 - 1) * uvec/(gamma(m/2 + 1) * 2 ^ 
                  ((m - 1)*(m - 2)/2 + 1) * prod(bvec)))^(1/m)
  lng <- sqrt(apply(tmp, 1, crossprod))
  out <- tmp*rvec/lng
  out
}