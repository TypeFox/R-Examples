#' Number of nodes at each depth
#' 
#' Find the number of nodes at each depth in the tree
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @return A vector of widths, where entry i is the number of nodes at depth i. There is a single node at depth 0 (the root) which is not included in the vector, for simplicity.
#' 
#' @import ape
#'   
#' @examples
#' ## Find the node widths in a random tree with 10 tips:
#' tree <- rtree(10)
#' tree$edge.length <- rep(1,18) # to make it easier to see the width and depths in the plot
#' plot(tree)
#' widths(tree)
#' 
#' @export
widths <- function(tree) {
  depths <- getDepths(tree)
  allDepths <- c(depths$tipDepths,depths$nodeDepths)
  widths <- sapply(1:max(allDepths), function(x) sum(allDepths==x))
  names(widths) <- 1:max(allDepths)
  return(widths)
}