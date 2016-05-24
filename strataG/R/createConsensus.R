#' @title Consensus Sequence
#' @description Return a consensus sequence from set of aligned sequences, 
#'   introducing IUPAC ambiguity codes where necessary.
#' 
#' @param x a \linkS4class{gtypes} object with aligned sequences or a list of 
#'   aligned DNA sequences.
#' @param ignore.gaps logical. Ignore gaps at a site when creating consensus. 
#'   If \code{TRUE}, then bases with a gap are removed before consensus is 
#'   calculated. If \code{FALSE} and a gap is present, then the result is a gap.
#'   
#' @return A character vector of the consensus sequence.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.seqs)
#' createConsensus(dolph.seqs)
#' 
#' @export
#' 
createConsensus <- function(x, ignore.gaps = FALSE) { 
  x <- as.multidna(x)
  
  result <- lapply(getSequences(x, simplify = FALSE), function(dna) {
    dna <- as.character(as.matrix(dna))
    apply(dna, 2, iupacCode, ignore.gaps = ignore.gaps)
  })
  
  if(length(result) == 1) {
    result[[1]]
  } else {
    names(result) <- getLocusNames(x)
    result
  }
}