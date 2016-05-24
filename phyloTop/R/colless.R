#' Colless number
#' 
#' Finds the Colless number for a tree.
#' Note that the package \code{apTreeshape} has a function \code{colless} to compute the Colless imbalance with additional options to normalise it based on the model;
#' we include this simple function here for convenience within this package, and for use on objects of class \code{phylo} and \code{phylo4}.
#' 
#' @author Michael Boyd \email{mboyd855@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param normalise option to normalise the result by dividing by the number of tip pairs. Defaults to \code{TRUE}.
#' @return The Colless imbalance number of the tree.
#' 
#' @import ape
#'   
#' @examples
#' ## Find the Colless imbalance of a random tree with 10 tips:
#' tree <- rtree(10)
#' plot(tree)
#' colless.phylo(tree)
#' 
#' 
#' @export
colless.phylo <- function(tree,normalise=TRUE) {
  ntips <- length(tree$tip.label)
  if (ntips==2) {return(0)}
  tImb <- treeImb(tree)
  diffs <- abs(apply(tImb,1,diff))
  if (normalise) {
    n <- ((ntips-1)*(ntips-2))/2
    return(sum(diffs)/n)
  }
  return(sum(diffs))
}
