#' Evaluate integral
#'
#' This function will evaluate the integral in equation 
#' (4) on page 7 of the paper. See Details section.
#' 
#' @param limit.vec This is a vector of length 2. It consists of the lower and
#' upper limits in the integral. Checking is carried out to ensure that is of
#' length two, and that limit.vec[1] <= limit.vec[2].
#' @param theta.diff A vector of length p-1, where p is the number of populations
#' of treatments. Coordinate [i] in theta.diff corresponds to \eqn{\theta_i -
#' \theta_{i+1}}. See \code{\link{genDelMat}}.
#' @param sigma.2 The known variance of the error terms.
#' @param n The number of replications per population.
#'
#' @export
#'
#' @details This function evaluates the integral, and works with the lower and
#' upper limits that it is given. If one desires to compute the coverage
#' probability for an interval defined by \eqn{X_{(1)} \pm c}, then the user should
#' look at the function \link{exactCoverageProb} in this package.
#'
#' @examples
#' del1 <- c(2, 4)
#' integrate2(c(-1.1,1.3), del1)
#'
#' @seealso
#' \link{exactCoverageProb}, \link{integrand}
#'
#' @return The function returns a scalar value that is the value of the integral
#' in equation (4) of page 7, defined by the lower and upper limits provided
#' here.

integrate2 <- function(limit.vec=c(-3,3), theta.diff, sigma.2 = 1, n=1) {
  # error checking:
  if(length(limit.vec)!=2)
    stop("limit.vec should be a vector of length 2")
  if(limit.vec[2] < limit.vec[1])
    stop("limit.vec must have limit.vec[1] <= limit.vec[2]")

  mat1 <- genDelMat(theta.diff, sigma.2, n)
  diag(mat1) <- NA

  integrate(integrand, lower=limit.vec[1], upper=limit.vec[2], mat1=mat1)$value
}
