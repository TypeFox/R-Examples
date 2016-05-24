#' Array data of Generalized Procrustes Analysis (GPA) coordinates of dorsal anchors of \emph{Ligophorus liewi} samples
#'
#' An array data of GPA coordinates of dorsal anchors (left-right averaged) from 31 \emph{Ligophorus liewi} samples.
#' @docType data
#' @usage data(liewi_gpa)
#' @format An array of 31 matrices, each having 11 rows (landmarks) and 2 columns (xy GPA coordinates)
#' @details All samples in this array have quality score above 10.
#' @keywords datasets
#' @source Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' @seealso \code{\link{Qscore}}, \code{\link{plotLM}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @examples
#' data(liewi_gpa)
#'
#' nice_title <- expression(paste("Dorsal anchor ", italic(L.liewi)))
#' plotLM(liewi_gpa, tit=nice_title, pointscale=0.8, axispointscale=0.8, 
#' meansize=1.2, polygon.outline=TRUE,c(-.6,.6),c(-.6,.6) )
#'
"liewi_gpa"
