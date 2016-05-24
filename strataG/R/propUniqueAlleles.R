#' @title Proportion Unique Alleles
#' @description Calculate the proportion of alleles that are unique.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a vector of the proportion of unique (occuring only in one individual) 
#'   alleles for each locus.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.seqs)
#' strata <- dolph.strata$fine
#' names(strata) <- dolph.strata$ids
#' dloop <- sequence2gtypes(dolph.seqs, strata, seq.names = "dLoop")
#' dloop <- labelHaplotypes(dloop)$gtypes
#' 
#' propUniqueAlleles(dloop)
#' 
#' @export
#' 
propUniqueAlleles <- function(g) { 
  id.rows <- sub("\\.[[:digit:]]*$", "", rownames(loci(g)))
  apply(loci(g), 2, function(locus) {
    id.a.freqs <- table(id.rows, locus)
    is.unique <- apply(id.a.freqs, 2, function(x) sum(x > 0) == 1)
    sum(is.unique)
  }) / numAlleles(g)
}