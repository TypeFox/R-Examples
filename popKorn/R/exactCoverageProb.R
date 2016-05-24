#' Evaluate exact coverage probability
#'
#' This function will evaluate the exact coverage probability, as given in
#' equation (4) on page 7 of the paper. See Details section.
#' 
#' @param c.vec This is a vector of length 2. It consists of the lower and
#' upper limits in the integral. Checking is carried out to ensure that is of
#' length two, and that 0 <= c.vec[1] <= c.vec[2]. This parameter is ignored if
#' \code{lambda} is not missing.
#' @param theta.diff A vector of length p-1, where p is the number of populations
#' of treatments. Coordinate [i] in theta.diff corresponds to \eqn{\theta_i -
#' \theta_{i+1}}. See \code{\link{genDelMat}}.
#' @param lambda In case the user wishes to use the shrinkage version, this
#' parameter should be specified. It must be between 0 and 1. 
#' @param c.val In case lambda is specified, this must not be missing. This will
#' be combined with lambda to create a c.vec. This very function will then call
#' itself.
#' @param sigma.2 The known variance of the error terms.
#' @param n The number of replications per population.
#'
#' @export
#'
#' @details This function evaluates the coverage probability for an interval defined 
#' by \eqn{(X_{(1)} - c_2, X_{(1)} + c_1)}. Note that, as specified in the
#' reference paper, we must have that \eqn{0 \le c_1 \le c_2}. This function
#' will call \link{integrate2}. Please note the ordering of the elements in the
#' \code{c.vec} argument: the first element corresponds to the upper limit of
#' the interval, and to the negative of the lower limit of the integral.
#'
#' @examples
#' del1 <- c(2, 4)
#' exactCoverageProb(c(1.1,1.3), del1)
#' exactCoverageProb(theta.diff=c(2,3,4), lambda=0.9, c.val=2)
#'
#' @seealso
#' \link{integrate2}, \link{integrand}
#'
#' @return The function returns a scalar value that is the value of the exact
#' coverage coverage probability defined in equation (4) of page 7.

exactCoverageProb <- function(c.vec, theta.diff, lambda, c.val, sigma.2=1,
                              n=1) {
  # error checking:
  if(!missing(lambda)) {
    if(lambda < 0 || lambda > 1)
      stop("lambda parameter must be between 0 and 1.")
    if(missing(c.val))
      stop("If lambda is specified, c.val must also be specified.")
    if(c.val <= 0)
      stop("c.val must be positive")
    c.vec.rec <- c(lambda*c.val, c.val)
    return(exactCoverageProb(c.vec.rec, theta.diff))
  }
  if(length(c.vec)!=2)
    stop("c.vec should be a vector of length 2")
  if(any(c.vec < 0))
    stop("c.vec should be non-negative")
  if(c.vec[2] < c.vec[1])
    stop("c.vec must have c.vec[1] <= c.vec[2]")

  integrate2(limit.vec=c(-c.vec[1], c.vec[2]), theta.diff, sigma.2, n)
}
