#' @title Nucleotide Diversity
#' @description Calculate nucleotide diversity for set of haplotypes.
#' 
#' @param x a set of sequences or a \linkS4class{gtypes} object with sequences.
#' @param bases nucleotides to consider when calculating diversity.
#' 
#' @return Nucleotide diversity by site. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.seqs)
#' strata <- dolph.strata$fine
#' names(strata) <- dolph.strata$ids
#' dloop <- sequence2gtypes(dolph.seqs, strata, seq.names = "dLoop")
#' 
#' nucleotideDiversity(dloop)
#' 
#' @importFrom swfscMisc diversity
#' @export
#' 
nucleotideDiversity <- function(x, bases = c("a", "c", "g", "t")) {
  x <- as.multidna(x)
  bases <- tolower(bases)
  result <- lapply(getSequences(x, simplify = FALSE), function(dna) {  
    dna <- as.character(as.matrix(dna))
    site.div <- apply(dna, 2, function(b) {
      swfscMisc::diversity(b[b %in% bases])
    })
    names(site.div) <- 1:length(site.div)
    site.div
  })
  
  if(length(result) == 1) {
    result[[1]]
  } else {
    names(result) <- getLocusNames(x)
    result
  }
}