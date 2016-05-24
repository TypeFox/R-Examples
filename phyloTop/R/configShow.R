#' Plot a tree highlighting configurations
#' 
#' Plot a tree, highlighting configurations of a given size.
#' 
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#' @author Michael Boyd \email{mboyd855@gmail.com}
#'   
#' @param tree a tree of class \code{phylo4}
#' @param configSize an integer giving the configuration size of interest
#' @param mainCol colour for branches which are not in configurations of the chosen size (default is black)
#' @param configCol colour for branches which are in such configurations (default is red)
#' @param ... further arguments to be passed to plot.phylo
#' @return A plot of the tree, highlighting the configurations of the given size.
#' 
#' @import ape
#' 
#' @examples
#' ## Highlight pitchforks in a random tree with 20 tips:
#' configShow(rtree(20),3, edge.width=2)
#' 
#' @export
configShow <- function(tree, configSize, mainCol="black", configCol="red",  ...) {
  tree <- phyloCheck(tree)
  edgeList <- tree$edge
  nEdges <- length(tree$edge[,1])
  col <- rep(mainCol,nEdges)
  allcladeSizes <- nConfig(tree)$cladeSizes
  cladesToHighlight <- which(allcladeSizes==configSize)
  
  ntips=length(tree$tip.label)
  nn=tree$Nnode
  Ancs=(ntips+1):(ntips+nn) # assumes tips are 1:ntip, then internal nodes
  
  # for each internal node (numbered 1:(ntips-1) in the rows of "pointers"), find its immediate children 
  Pointers=t(vapply(Ancs, function(x) tree$edge[tree$edge[,1]==x,2], FUN.VALUE=c(1,2))) 
  
  for (i in cladesToHighlight){
    col[which(edgeList[,1]==i)] <- configCol
    node <- i
        while(any(is.element(node, Ancs))) {# while any of the nodes are ancestors (i.e. have descendants)
          children <- Pointers[(node[which(is.element(node,Ancs))] - ntips),] 
          col[which(edgeList[,1]%in%children)] <- configCol
          node <- children
        }
  }

  plot.phylo(tree, edge.color=col, ...)
}