#' Node depth
#' 
#' Determine the depth of a particular node in a tree, defined as the number of edges between it and the root.
#' (So the root has depth zero, its children have depth one, etc.)
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param node a number corresponding to a node in the tree.
#' 
#' @return An integer corresponding to the depth of the given node.
#' 
#' @import ape
#'   
#' @examples
#' ## Find the depth of node 34 in a random tree with 20 tips:
#' tree <- rtree(20)
#' plot(tree)
#' nodelabels()
#' nodeDepth(tree,34)
#'  
#' 
#' @export
nodeDepth <- function(tree,node) {
  depths <- getDepths(tree)
  allDepths <- c(depths$tipDepths,depths$nodeDepths)
  l <- length(allDepths)
  if (node > l) stop(paste0("Please supply a valid node number between 1 and ",l))
  return(allDepths[[node]])
}
