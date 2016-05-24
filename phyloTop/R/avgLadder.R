#' Average ladder size
#'
#' Finds the mean size of ladders in the tree
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param normalise option to normalise the result, default is \code{FALSE}
#' @return The mean ladder size
#' 
#' @import ape
#' 
#' @seealso \code{\link{ladderSizes}}
#'
#' @examples
#' ## Find the average ladder size in a random tree with 20 tips:
#' tree <- rtree(20)
#' plot(tree)
#' avgLadder(tree)
#' # and the normalised average ladder size:
#' avgLadder(tree, normalise=TRUE)
#' 
#'
#' @export
avgLadder <- function(tree, normalise=FALSE) {
  l <- ladderSizes(tree)$ladderSizes
  if (normalise==FALSE) {return(mean(l))}
  else {
    ntips <- length(tree$tip.label)
    if (ntips==2) {return(0)}
    else return(mean(l)/(ntips-2))
    }
}