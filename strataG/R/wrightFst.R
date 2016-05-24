#' @title Wright's Fst
#' @description Calcualte Wright's Fst from Ne, dispersal, and generation time.
#' 
#' @param Ne Effective population size.
#' @param dispersal migration rate in terms of probability of an individual 
#'   migrating in a generation.
#' @param gen.time number of generations since ancestral population.
#' @param ploidy ploidy of the locus
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
wrightFst <- function(Ne, dispersal, gen.time, ploidy) {
  1 / (2 * ploidy * Ne * dispersal * gen.time + 1)
}