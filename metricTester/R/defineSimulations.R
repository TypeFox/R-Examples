#' Output all spatial simulations as a list of named functions
#'
#' Creates a list of named functions, each of which accept a simulations.input object
#'
#' @details All of the spatial simulations we calculated for our manuscript are included 
#' in this function. To add additional spatial simulations, it can either be defined on
#' the fly or to permanently include a new simulation in all downstream simulations, it
#' can be included here. The function needs to be included with a name, and it must accept
#' a simulations.input object. If the function needs additional elements not included in 
#' that input, then the prepSimulations function must also be revised.
#'
#' @return A list of named functions
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' defineSimulations()

defineSimulations <- function()
{
	list("random"=randomArena, "filtering"=filteringArena, "competition"=competitionArena)
}
