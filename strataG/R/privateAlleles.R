#' @title Private Alleles
#' @description The number of private alleles in each strata and locus.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return matrix with the number of private alleles in each strata at each locus.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' privateAlleles(msats.g)
#' 
#' @export
#' 
privateAlleles <- function(g) {
  freqs <- alleleFreqs(g, TRUE)
  
  do.call(rbind, lapply(freqs, function(f) {
    f <- f[, "freq", , drop = FALSE]
    f[f > 0] <- 1
    pa <- rowSums(apply(f, 1, function(x) {
      if(sum(x > 0) == 1) x else rep(0, length(x))
    }))
    names(pa) <- dimnames(f)[[3]]
    pa
  }))
}