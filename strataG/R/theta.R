#' @title Theta
#' @description Calculate theta from heterozygosity of each locus.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return vector of theta values for each locus.
#' 
#' @details Calculates theta for each locus using the 
#'   \code{\link[pegas]{theta.h}} function.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom pegas theta.h
#' @importFrom stats na.omit
#' @export
#' 
theta <- function(g) {
  sapply(colnames(g@loci), function(x) theta.h(na.omit(g@loci[, x])))
}