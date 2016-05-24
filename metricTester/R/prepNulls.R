#' Prep data for null randomizations
#'
#' Given a phylo object, a picante-style community data matrix (sites are rows,
#' species are columns), and an optional vector of regional abundance, prepare data for
#' randomizations.
#'
#' @param tree Phylo object
#' @param picante.cdm A picante-style community data matrix with sites as rows, and
#' species as columns
#' @param regional.abundance A character vector in the form "s1, s1, s1, s2, s2, s3, etc".
#' Optional, will be generated from the input CDM if not provided.
#' @param distances.among An optional symmetric distance matrix describing the distances
#' among plots/etc, for use with null models like the dispersal null.
#' 
#' @details Returns a named list with four elements: the original phylogenetic tree,
#' the original picante-style CDM, a spacodi-style CDM, and vector of regional abundance.
#'
#' @return A list of class nulls.input
#'
#' @export
#'
#' @importFrom spacodiR as.spacodi
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
#' prepped <- prepNulls(tree, cdm)

prepNulls <- function(tree, picante.cdm, regional.abundance=NULL, distances.among=NULL)
{
	spacodi.cdm <- suppressMessages(as.spacodi(picante.cdm))
	if(is.null(regional.abundance))
	{
		warning("Regional abundance not provided. Assumed to be equivalent to CDM")
		regional.abundance <- abundanceVector(picante.cdm)
	}
	if(is.null(distances.among))
	{
		warning("Distances among plots not provided. Null models that require this input will not be run")
		distances.among <- "ignore"
	}
	dat <- list("tree"=tree, "picante.cdm"=picante.cdm, "spacodi.cdm"=spacodi.cdm, 
		"regional.abundance"=regional.abundance, "distances.among"=distances.among)
	class(dat) <- c("list", "nulls.input")
	dat
}
