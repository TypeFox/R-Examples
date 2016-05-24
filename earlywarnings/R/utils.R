#' @import Kendall
#' @import KernSmooth
#' @import lmtest
#' @import nortest
#' @import som
#' @import spam
#' @import stats
NULL

#' Description: Get group assigment indices for univariate data points, given cluster break points
#'
#' Arguments:
#'    @param x Univariate data vector
#'    @param breakpoints Cluster breakpoints
#'
#' Returns:
#'   @return A vector of cluster indices
#'
#' @export
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples #
#'
#' @keywords early-warning

UnivariateGrouping <- function(x, breakpoints) {
    g <- rep.int(NA, length(x))
    mps <- c(breakpoints, Inf)
    for (i in 1:length(mps)) {
        g[x <= mps[[i]] & is.na(g)] <- i
    }
    g
} 
