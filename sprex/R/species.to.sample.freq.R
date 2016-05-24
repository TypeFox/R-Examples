#' @title Create Vector of Sample Frequencies
#' @description Create vector of sample frequencies from vector of 
#'   species frequencies.
#' 
#' @param f a vector of species frequencies where \code{f[i]} is the 
#'   number of species represented by only \code{i} samples.
#' 
#' @return a vector(\code{x}) where \code{x[i]} is of the number of samples 
#'   in the \code{i}-th species.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @seealso sample.to.species.freq
#' 
#' @examples
#' data(osa.old.growth)
#' f <- expand.freqs(osa.old.growth)
#' x <- species.to.sample.freq(f)
#' print(f)
#' print(x)
#' 
#' @export

species.to.sample.freq <- function(f) {
  unlist(lapply(1:length(f), function(i) rep(i, f[i])))
}
