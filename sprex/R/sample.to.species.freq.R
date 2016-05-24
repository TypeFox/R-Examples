#' @title Create Vector of Species Frequencies
#' @description Create vector of species frequencies from vector of sample 
#'   frequencies.
#' 
#' @param x a vector where \code{x[i]} is of the number of samples in the 
#'   \code{i}-th species.
#' @param min.f minimum size of return vector. Return vector is zero-padded up 
#'   to this length if it would normally be shorter.
#' 
#' @return a vector(\code{f}) of species frequencies where \code{f[i]} is the 
#'   number of species represented by only \code{i} samples.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @seealso species.to.sample.freq
#' 
#' @examples
#' x <- sample(1:20, 20, rep = TRUE)
#' f <- sample.to.species.freq(x)
#' print(x)
#' print(f)
#' 
#' @export
#' 
sample.to.species.freq <- function(x, min.f = NULL) {
  x.df <- as.data.frame(table(x), stringsAsFactors = FALSE)
  x.df[, 1] <- as.numeric(as.character(x.df[, 1]))
  x.df <- x.df[x.df[, 1] > 0, ]
  rownames(x.df) <- NULL
  f <- expand.freqs(x.df)
  if(!is.null(min.f)) c(f, rep(0, max(0, min.f - length(f)))) else f
}
