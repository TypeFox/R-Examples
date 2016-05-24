#' Show ladders
#' 
#' Plot a tree, highlighting any ladders within it.
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#' @author Michael Boyd \email{mboyd855@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param mainCol colour for edges which are not ladders (default is black)
#' @param ladderEdgeCol colour for ladder edges (default is red)
#' @param ladderNodeCol colour for ladder nodes (default is red)
#' @param ... further arguments to be passed to plot.phylo
#' 
#' @return A plot of the tree, with ladder edges and nodes highlighted by colour.
#' 
#' @seealso \code{\link{ladderSizes}}
#'   
#' @import ape
#' 
#' @examples
#' ## Highlight in blue the ladders in a random tree with 50 tips:
#' tree <- rtree(50)
#' ladderShow(tree, edge.width=2)
#' # compare to:
#' ladderSizes(tree)
#' 
#' 
#' @export
ladderShow <- function(tree, mainCol="black", ladderEdgeCol="red", ladderNodeCol="red", ...) {
  tree <- phyloCheck(tree)
  edgeList <- tree$edge
  nEdges <- nrow(tree$edge)
  col <- rep(mainCol,nEdges)
  ladderSizes <- ladderSizes(tree)
  
  ladderEdges <- ladderSizes$ladderEdges
  ladderNodes <- (1:(2*length(tree$tip.label)-1))[ladderSizes$ladderNodes] # coerces into correct "integer" format

  # If the ladder number is at least one, change the edge colour
  for (i in ladderEdges) {
    col[i] <- ladderEdgeCol
  }
  plot.phylo(tree,edge.color=col, ...)
  if(length(ladderNodes)!=0) nodelabels(node=ladderNodes, bg=ladderNodeCol)
}


