#' Prep data to test phylogenetic community structure metrics
#'
#' Given a phylo object, and a picante-style community data matrix (sites are rows,
#' species are columns), prepare data for analysis.
#'
#' @param tree Phylo object. 
#' @param picante.cdm A picante-style community data matrix with sites as rows, and
#' species as columns
#' @param optional.dists A symmetric distance matrix can be directly supplied. This option
#' is experimental. Behavior depends on metric being used. If the metric in question
#' relies on the phylogenetic distance matrix from a call to cophenetic(tree), then this 
#' optional distance matrix will be inserted instead. 
#' 
#' @details Returns a named list with three elements: the original phylogenetic tree
#' phylogenetic distances among species, and the original picante-style CDM.
#'
#' @return An object of class metrics.input
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
#'
#' prepped <- prepData(tree, cdm)

prepData <- function(tree, picante.cdm, optional.dists=NULL)
{
	if(!is.null(optional.dists))
	{
		#this is a very rudimentary check to see if the optional distance matrix contains
		#all spp that are in the CDM. note that it can contain species that are not in the
		#CDM, but not vice versa
		if(length(setdiff(colnames(picante.cdm), row.names(optional.dists))) > 0)
		{
			stop("Some species in your CDM are not in your distance matrix")
		}
		dists <- optional.dists
	}
	else
	{
		dists <- ape::cophenetic.phylo(tree)
	}
		
	#you have a check in plotContents to exclude plots w < 2 spp, but after
	#randomizations it is possible to end up with plots that include < 2 spp. exclude
	#these. note that dplyr does not need even sample sizes or anything like that, so
	#this should hopefully work
	picante.cdm <- picante.cdm[apply(picante.cdm, 1, lengthNonZeros) >= 2,]
	dat   <- list("tree"=tree, "dists"=dists, "picante.cdm"=picante.cdm)
	class(dat) <- c("list", "metrics.input")
	dat
}
