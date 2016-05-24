#' @name contCA00
#' @title Simulated 7 state (1+1)D field
#' @docType data
#' @usage 
#' data(contCA00)
#' @format
#' Contains the running example dataset used in hard LICORS & mixed LICORS.
#' 
#' A list with three \eqn{(1+1)D} fields, each one extending over \eqn{N = 100} 
#' pixels in space, and \eqn{T = 200} over time:
#' \itemize{
#'   \item \code{observed}
#'   \item \code{states}
#'   \item \code{predictive_states}
#' }
#' 
#' @keywords dataset
#' @references \url{arxiv.org/abs/1206.2398}
#' @examples
#' # set original par parameters
#' op <- par(no.readonly=TRUE)
#' 
#' data(contCA00)
#' par(mfrow = c(2,2), mar = c(3,3,2,1))
#' for (ii in 1:3){
#'   image2(contCA00[[ii]], legend = FALSE, col ="RdBu", 
#'          main = attr(summary(contCA00), "dimnames")[[1]][ii])
#'   mtext("Time", 1, 1)
#'   mtext("Space", 2, 1)
#' }
#' par(op)
#' \dontrun{
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 2, FLC = 0), shape ="cone")
#' bb = data2LCs(contCA00$observed, LC.coordinates = LC_geom$coordinates)
#' image2(bb$PLC)
#' image2(cor(bb$PLC), zlim = c(-1,1), col = "RdBu")
#' mod_kk = kmeanspp(bb$PLC, k = 10)
#' plot(bb$FLC, col = mod_kk$cluster, pch = ".", cex = 3)
#' 
#' ff = estimate_LC_pdfs(bb$FLC, states = mod_kk$cluster, method = "nonparametric")
#' matplot(bb$FLC, ff, pch = ".", cex = 2)
#' }
#' 
NULL
