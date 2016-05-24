#' @title Permutation Test P-value
#' @description Calculate the p-value for a permutation test.
#'
#' @param obs observed value.
#' @param null.dist vector of values from permutation null distribution.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' null.dist <- rnorm(1000)
#' obs <- rnorm(1, mean = 1)
#'
#' plot(density(null.dist), xlim = range(c(obs, null.dist)), main = "")
#' abline(v = obs)
#' print(obs)
#' pVal(obs, null.dist)
#'
#' @export
#' 
pVal <- function(obs, null.dist) {
  (sum(null.dist >= obs) + 1) / (length(null.dist) + 1)
}