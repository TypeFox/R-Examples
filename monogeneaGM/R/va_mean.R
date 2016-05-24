#' Mean Generalized Procrustes Analysis (GPA) coordinates for ventral anchors of 13 \emph{Ligophorus} species
#'
#' A data set containing mean GPA coordinates for the ventral anchors of 13 \emph{Ligophorus} species.
#' @docType data
#' @usage data(va_mean)
#' @format A list containing 13 objects, each object being a matrix of 11 rows (landmarks) and
#' 2 columns (xy-coordinates) representing GPA-coordinates averaged from samples belonging to
#' the same species
#' @details For each specimen, GPA coordinates of all landmarks are averages of those from 
#' left and right anchors. Mean GPA coordinates of landmarks are required for shape evolution analysis as well as
#' performing Adam's Kmult test for presence of phylogenetic signal (Adams, 2014).
#' @keywords datasets
#' @source Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @seealso \code{\link{shapeEvo}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Adams DC. (2014). A generalized K statistic for estimating phylogenetic signal
#' from shape and other high-dimensional multivariate data. Systematic Biology 63: 685-697.
"va_mean"