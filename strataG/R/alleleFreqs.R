#' @title Allele Frequencies
#' @description Calculate allele frequencies for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical. If \code{TRUE} every element in the return list is 
#'   a three dimensional array where the third dimension contains frequencies 
#'   and proportions for each stratum.
#'
#' @return A list of allele frequencies for each locus. Each element is a
#'   matrix or array with frequencies by count (\code{freq}) and 
#'   proportion (\code{prop}) of each allele.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' f <- alleleFreqs(msats.g)
#' f[[1]]
#' 
#' f.pop <- alleleFreqs(msats.g, TRUE)
#' f.pop[[1]]
#' 
#' @export
#' 
alleleFreqs <- function(g, by.strata = FALSE) {
  freqs <- vector("list", length = ncol(g@loci))
  strata <- rep(g@strata, g@ploidy)
  for(i in 1:length(freqs)) {
    if(by.strata & nlevels(strata) != 1) {
      f <- table(g@loci[, i], strata)
      p <- prop.table(f, 2)
      freqs[[i]] <- array(dim = list(nrow(f), 2, ncol(f)))
      for(j in 1:ncol(f)) freqs[[i]][, , j] <- cbind(f[, j], p[, j])
      dimnames(freqs[[i]]) <- list(rownames(f), c("freq", "prop"), colnames(f))
    } else {
      f <- table(g@loci[, i])
      freqs[[i]] <- cbind(freq = f, prop = f / sum(f))
    }
  }
  names(freqs) <- colnames(g@loci)
  freqs
}
