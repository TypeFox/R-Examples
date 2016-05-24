#' @title Plot the species abundance distribution (SAR), i.e. objects of class sar
#'
#' @description Plot species or endemics area relationship with flexibility to adjust plotting parameters
#'
#' @details see examples
#' 
#' 
#' @param x an object of class SAR made with 
#' @param add logical; should new \code{sar} object be added to current plot or made its own plot
#' @param ... arguments passed to \code{plot}
#' 
#' @export
#' 
#' @importFrom graphics plot points
#' @examples
#' data(anbo)
#' anbo.obs.sar <- empiricalSAR(anbo$spp, anbo$count, anbo$row, anbo$col, Amin=1, A0=16)
#' plot(anbo.obs.sar)

#' @return NULL
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso empiricalSAR, downscaleSAR, upscaleSAR, meteSAR
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

plot.sar <- function(x, add=FALSE, ...) {
    if(add) {
    	points(x[['A']], x[['S']], type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
    } else {
        plot(x[['A']], x[['S']], xlab='Area', 
             ylab=sprintf('Number of %s', ifelse(attr(x, 'type') == 'ear', 'endemics', 'species')), 
             type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
    }
}
