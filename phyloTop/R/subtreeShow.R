#' Highlight a subtree
#' 
#' Plot a tree, highlighting the clade(s) descending from the given node(s)
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @param nodeList a list of one or more internal nodes in the tree.
#' @param showNodeLabels option of whether to show node labels. Default is "nodeList", which only labels the nodes in \code{nodeList}.
#' Choosing \code{showNodeLabels="all"} will display all node labels; any other arguments will remove all node labels.
#' @param mainCol colour for the edges which are not highlighted (default is black).
#' @param subtreeCol colour for the edges which are highlighted (default is red).
#' @param nodeLabelCol background colour for any node labels shown (default is light blue)
#' @param ... further arguments to be passed to plot.phylo
#' 
#' @return A plot of the tree, with clade(s) descending from the given node(s) highlighted.
#' 
#' @import ape
#'   
#' @examples
#' ## Highlight the clade(s) descending from nodes 23 and 35 in a random tree on 20 tips:
#' tree <- rtree(20)
#' subtreeShow(tree, nodeList=c(23,35))
#' # change aesthetics:
#' subtreeShow(tree,nodeList=c(23,35), mainCol="navy", subtreeCol="cyan", 
#'    nodeLabelCol="cyan", edge.width=2)
#'    
#' 
#' @export
subtreeShow <- function(tree, nodeList, showNodeLabels="nodeList", mainCol="black", subtreeCol="red", nodeLabelCol="lightblue", ...) {
  tree <- phyloCheck(tree)
  ntips <- length(tree$tip.label)
  if (any(nodeList %in% 1:ntips)) stop("Each node in nodeList needs to be an internal node, numbered in the range ",(ntips+1),":",(2*ntips-1))
  
  originalNodeList <- nodeList # for later colouring
  edgeList <- tree$edge
  edgeCol <- rep(mainCol,length(edgeList))
  
  # for each node of nodeList, which two rows of edgeList are its immediate descendant branches
  edgesToColour <- which(edgeList[,1] %in% nodeList)
  
  # initialise "tmp" which keeps adding edges to be coloured in the following while loop
  tmp <- edgesToColour
  
  while (any(nodeList %in% (ntips+1):(2*ntips-1))) {
  # update nodeList to be the nodes at the ends of these edges
  nodeList <- edgeList[tmp,2]
  tmp <- which(edgeList[,1] %in% nodeList)
  edgesToColour <- c(edgesToColour,tmp)
  }
  
  edgeCol[edgesToColour] <- subtreeCol
  plot.phylo(tree, edge.color=edgeCol, ...)
  if (showNodeLabels=="nodeList") {nodelabels(node=originalNodeList, bg=nodeLabelCol)}
  else if (showNodeLabels=="all") {nodelabels(bg=nodeLabelCol)}
}

