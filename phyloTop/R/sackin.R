#' Sackin index
#' 
#' Finds the Sackin index for a tree. 
#' Note that the package \code{apTreeshape} has a function \code{sackin} to compute the Sackin index with additional options to normalise it based on the model;
#' we include this simple function here for convenience within this package, and for use on objects of class \code{phylo} and \code{phylo4}.
#' 
#' @author Michael Boyd \email{mboyd855@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param normalise option to normalise the result, default is \code{FALSE}.
#' 
#' @return The Sackin index of the tree.
#' 
#' @import ape
#'   
#' @examples
#' ## Sackin index of a random tree with 10 tips:
#' sackin.phylo(rtree(10))
#' 
#' ## normalised Sackin index:
#' sackin.phylo(rtree(10), normalise=TRUE)
#' 
#' 
#' @export
sackin.phylo <- function(tree,normalise=FALSE) {
  depths <- getDepths(tree)
  tipDepths <- depths$tipDepths - 1
  if (normalise==FALSE){
  return(sum(tipDepths))
  }
  else {
  n <- length(tree$tip.label)
  
  return((sum(tipDepths))/((1/2)*(n*(n+1)) - 1)) 
  }
}