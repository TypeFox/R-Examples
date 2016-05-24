#' Number of pitchforks
#' 
#' Finds the number of pitchforks in a tree. A pitchfork is considered to be a clade with three tips.
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'     
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param normalise option to normalise the result, default is \code{FALSE}.
#' 
#' @return An integer representing the number of pitchforks in the tree.
#'  
#' @seealso \code{\link{configShow}}   
#'     
#' @import ape
#' 
#' @examples
#' ## Find the number of pitchforks in a random tree with 20 tips:
#' tree <- rtree(20)
#' plot(tree)
#' pitchforks(tree)
#' # and the normalised pitchfork number:
#' pitchforks(tree, normalise=TRUE)
#'
#' ## Note that the function configShow can be used to highlight the pitchforks in the tree:
#' configShow(tree, 3, edge.width=2)
#' 
#' 
#' @export
pitchforks<-function(tree, normalise=FALSE) {
  tree <- phyloCheck(tree)
  ntips <- length(tree$tip.label)
  if (ntips==2) {return(0)}
  if (normalise==FALSE) {return(nConfig(tree)$numClades[[3]])}
  else {return(3*nConfig(tree)$numClades[[3]]/ntips)}
}
