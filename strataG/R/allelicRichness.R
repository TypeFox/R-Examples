#' @title Allelic Richness
#' @description Calculate allelic richness for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return the allelic richness of each locus calculated as the number of 
#'   alleles divided by the number of samples without missing data at 
#'   that locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples 
#' data(msats.g)
#' allelicRichness(msats.g)
#'
#' @importFrom stats na.omit
#' @export
#' 
allelicRichness <- function(g) {
  apply(g@loci, 2, function(locus) {
    locus <- na.omit(locus)
    length(unique(locus)) / (length(locus) / ploidy(g))
  })
}
