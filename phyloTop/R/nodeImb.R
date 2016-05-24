#' Node imbalance
#' 
#' For a given node, this function gives the number of tips descending from each of its two children, as a measure of imbalance.
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param node a node index between 1 and 2n-1, where n is the number of tips in the tree.
#' 
#' @return Two integers corresponding to the number of tip descendants of each of the node's two children. If the node is itself a tip, then the vector (0,0) will be returned.
#' 
#' @seealso \code{\link{treeImb}}
#' 
#' @import ape
#'   
#' @examples
#' ## Find the imbalance of node 16 in a random tree with 10 tips:
#' tree <- rtree(10)
#' plot(tree)
#' nodelabels()
#' nodeImb(tree,16)
#' 
#' 
#' @export
nodeImb <- function(tree,node) {
  treeImb <- treeImb(tree)
  treeImb[node,]
}
