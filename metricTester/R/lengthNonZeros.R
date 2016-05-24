#' Calculate the species richness of a vector from a CDM
#'
#' Given a vector of abundances or presence/absences from a community data matrix, will
#' calculate the species richness of that vector.
#'
#' @param input.vector A vector from a community data matrix of abundances.
#' 
#' @details An internal function to calculate richness of a cdm. 
#'
#' @export
#'
#' @return A named vector of species richness. 
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #simulate tree with birth-death process
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' #note that with this example, each community in the cdm will be labeled by its richness
#' apply(cdm, 1, lengthNonZeros)

lengthNonZeros <- function(input.vector)
{
	nonZeros <- input.vector[input.vector != 0]
	return(length(nonZeros))
}
