#' @title Number of Generations to Equilibrium
#' @description Calculate the number of generations to equilibrium based on a an
#'   ideal Wright model.
#'
#' @param fst estimate of Fst between populations.
#' @param Ne population size.
#' @param gen.time generation time.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#' 
numGensEq <- function(fst, Ne, gen.time) {
  term1 <- log(1 - fst)
  term2 <- 1 - (1 / (2 * Ne))
  n.gens <- term1 / log(term2)
  eq.fst <- 1 - (term2 ^ gen.time)
  cbind(n.gens = n.gens, eq.fst = eq.fst)
}