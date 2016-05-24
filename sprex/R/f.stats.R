#' @title Frequency Vector Statistics
#' @description Number of observed species and samples in species frequency 
#'   vector.
#' 
#' @param f a vector of species frequencies where \code{f[i]} is the number 
#'   of species represented by only \code{i} samples.
#' 
#' @return a vector of the number of observed species (\code{s.obs}), 
#'   and the total number of samples (\code{n}).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @examples
#' data(osa.second.growth)
#' f <- expand.freqs(osa.second.growth)
#' f.stats(f)
#' 
#' @export
#' 
f.stats <- function(f) {
  s.obs <- sum(f)
  n <- sum(1:length(f) * f)
  c(s.obs = s.obs, n = n)
}