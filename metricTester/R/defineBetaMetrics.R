#' Output all beta metrics as a list of named functions
#'
#' Creates a list of named functions, each of which accept a metrics.input object
#'
#' @details All of the beta metrics we calculated for our manuscript are included in this
#' function. To add additional functions, they can either be defined on the fly or, to 
#' permanently include a new metric in all downstream simulations, it can be included
#' here. The function needs to be included with a name, and it must accept a metrics.input
#' as input. If the function needs additional elements not included in that input, then
#' the prepData function must also be revised.
#'
#' @return A list of named functions
#'
#' @export
#'
#' @importFrom picante mntd
#' @importFrom spacodiR spacodi.calc
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' defineBetaMetrics()

defineBetaMetrics <- function()
{
	list("richness"=my_betaRichness,
	"total_abundance"=my_totalAbundance,
	"Ist"=my_Ist,
	"Pst"=my_Pst,
	"Bst"=my_Bst,
	"PIst"=my_PIst,
	"meanNAW_MPD"=mean_naw_mpd,
	"meanInter_MPD"=mean_inter_mpd,
	"meanIntra_MPD"=mean_intra_mpd,
	"meanComplete_MPD"=mean_complete_mpd,
	"meanNAW_MNTD"=mean_naw_mntd,
	"meanAW_MNTD"=mean_aw_mntd
	)
}
