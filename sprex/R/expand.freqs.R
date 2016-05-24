#' @title Expand Frequency Matrix
#' @description Expand a matrix or data.frame of species frequencies to full 
#'   vector.
#' 
#' @param freq.mat a two column matrix or data.frame where the first column is 
#'   the number of samples, and the second column is the number of species 
#'   represented by with that many samples.
#' 
#' @return a vector(\code{f}) of species frequencies where each element 
#'   (\code{f[i]}) is the number of species represented by only \code{i} 
#'   samples.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(osa.old.growth)
#' f <- expand.freqs(osa.old.growth)
#' f
#' 
#' @export
#' 
expand.freqs <- function(freq.mat) {
  i.vec <- setdiff(1:max(freq.mat[, 1]), freq.mat[, 1])
  if(length(i.vec) > 0) {
    new.mat <- data.frame(i.vec, 0)
    colnames(new.mat) <- colnames(freq.mat)
    freq.mat <- rbind(freq.mat, new.mat)
  }
  freq.mat[order(freq.mat[, 1]), 2]
}
