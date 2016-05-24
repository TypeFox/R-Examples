#' Estimated Generalized Procrustes Analysis (GPA) coordinates of ventral anchors of root ancestor
#'
#' Estimated GPA landmark configuration of the ventral anchors of the root ancestor of 13 \emph{Ligophorus} species.
#' @docType data
#' @usage data(estimated_ancestral_va)
#' @format A matrix of 11 rows (landmarks) and 2 columns (xy-coordinates) representing the 
#' estimated GPA landmark configuration of the ventral anchors of the root ancestor (left and right anchors averaged). 
#' @details The root ancestor's GPA landmark configuration is unknown, but can be estimated using
#' similar data from extant species. The estimation is done using \code{fastAnc} function in the 
#' \code{phytools} package (Revell, 2012). Root ancestor mean GPA coordinates of anchor landmarks are required for 
#' shape evolution analysis.
#' @keywords datasets
#' @source Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @seealso \code{\link{shapeEvo}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Revell LJ. (2012). phytools: An R package for phylogenetic comparative biology (and other things). 
#' Methods in Ecology and Evolution 3:217-223.
#' @examples
#' library(gplots)
#' library(circular)
#'
#' data(va_mean) 
#' data(estimated_ancestral_va)
#' data(spcolmap)
#'
#' cladeII <- spcolmap$species[spcolmap$host %in% "L.subviridis"]
#' shapeEvo(va_mean, estimated_ancestral_va, col.lab="dodgerblue", 
#' clade=cladeII, exfac=2, tit="Ventral anchors")
#' 
#' #Some journals want the title to be left-adjusted, so set tit="" and then:
#' #title("a)", adj=0)
"estimated_ancestral_va"