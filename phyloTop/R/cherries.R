#' Cherry number
#' 
#' Finds the number of cherries in a tree. A cherry is considered to be a pair of sister tips.
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'     
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param normalise option to normalise the result, default is \code{FALSE}.
#' 
#' @return An integer representing the number of cherries in the tree.
#' 
#' @seealso \code{\link{configShow}}
#'
#' @import ape
#'        
#' @examples
#' ## Find the number of cherries in a random tree with 10 tips:
#' tree <- rtree(10)
#' plot(tree)
#' cherries(tree)
#' # and the normalised cherry number:
#' cherries(tree, normalise=TRUE)
#' 
#' ## Note that the function configShow can be used to highlight the cherries in the tree:
#' configShow(tree, 2, edge.width=2)
#' 
#' 
#' @export
cherries<-function(tree, normalise=FALSE) {
  tree <- phyloCheck(tree)
  if (normalise==FALSE) {return(nConfig(tree)$numClades[[2]])}
  else {return(2*nConfig(tree)$numClades[[2]]/length(tree$tip.label))}
}
