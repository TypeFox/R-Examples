#' Convert absolute abundance matrix to relative abundance
#'
#' Simple utility function to convert an absolute abundance matrix to a relative abundance
#' matrix.
#'
#' @param picante.cdm Picante-style community data matrix with
#' communities/plots/plots/etc as rows and species as columns
#' @param tree Optional phylo object
#' 
#' @details This function converts species' absolute abundances in a given community (a
#' row in the input CDM) into relative abundances by dividing observed abundances by the
#' maximum abundance in that row. If a tree is provided, the function confirms that the
#' CDM is indeed in the correct format, otherwise it assumes it is formatted correctly
#' and proceeds accordingly.
#'
#' @return A relative abundance matrix otherwise identical to the input CDM.
#'
#' @export
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

relativeCDM <- function(picante.cdm, tree=NULL)
{
	if(is.null(tree))
	{
		warning("assuming your cdm has species as columns, sites as rows")
	}
	else
	{
		if(length(setdiff(colnames(picante.cdm),tree$tip.label)) > 0)
		{
			stop("your cdm must have species as columns, sites as rows")
		}
	}
	
	newCDM <- t(apply(picante.cdm, 1, relative))
	newCDM
}

relative <- function(vect)
{
	result <- vect/max(vect)
	result
}
