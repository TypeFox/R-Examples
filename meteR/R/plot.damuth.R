#' @title Plot the relationship between abundance and metabolic rate, i.e. objects of class damuth
#'
#' @description Plot abundance-metabolic rate relationship with flexibility to adjust plotting parameters
#'
#' @details see examples
#' 
#' 
#' @param x an object of class damuth
#' @param add logical; should new \code{damuth} object be added to current plot or made its own plot
#' @param ... arguments passed to \code{plot}
#' 
#' @export
#' @importFrom graphics plot points

#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(arth$spp, arth$count, arth$mass^0.75)
#' ebar1 <- ebar(esf1)
#' plot(ebar1)

#' @return NULL
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso empiricalSAR, downscaleSAR, upscaleSAR, meteSAR
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

plot.damuth <- function(x, add=FALSE, ...) {
  if(add) {
    points(x[['n']], x[['e']], type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
  } else {
    plot(x[['n']], x[['e']], xlab='Abundance', 
         ylab='Mean metabolic rate', 
         type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
  }
}
