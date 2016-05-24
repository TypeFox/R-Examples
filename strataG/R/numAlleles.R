#' @title Number of Alleles
#' @description Return the number of alleles for each locus.
#'
#' @param g a \code{\link{gtypes}} object.
#'
#' @return vector of number of alleles per locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' numAlleles(msats.g)
#'
#' @importFrom stats na.omit
#' @export
#' 
numAlleles <- function(g) {
  apply(g@loci, 2, function(locus) length(unique(na.omit(locus))))
}