#' Maximum tree height
#' 
#' Find the maximum height of tips in the tree.
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param normalise option to normalise the result, default is \code{FALSE}.
#' 
#' @return An integer giving the maximum height of tips in the tree.
#' 
#' @import ape
#' 
#' @examples
#' ## Maximum height of tips in a random tree with 10 tips:
#' tree <- rtree(10)
#' plot(tree)
#' maxHeight(tree) 
#' maxHeight(tree, normalise=TRUE)
#' 
#' @export
maxHeight<-function(tree, normalise=FALSE){
  heights <- getDepths(tree)$tipDepths
  if (normalise==FALSE) {return(max(heights))}
  else {return(max(heights)/(length(tree$tip.label) - 1))}
}

