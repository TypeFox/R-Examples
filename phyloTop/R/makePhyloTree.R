#' Make phylogenetic tree (internal)
#' 
#' This function is called by \code{\link{getLabGenealogy}} to turn an edge list and corresponding lengths and root node number into a tree of class \code{phylo}
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param Edges a two-column matrix representing the edges of a tree. (It corresponds to the $edge entry of a \code{phylo} object, although the rows may not be in an expected order.) In each row the first entry is an internal node and the second entry is one of its immediate descendants.
#' @param Lengths a numeric vector of the same length as the columns of \code{Edges}. Each entry corresponds to the length of that edge.
#' @param Root the number of the internal node which is to be the root (usually the number of tips + 1)
#'
#' @return An object of class \code{phylo}. When called by \code{\link{getLabGenealogy}} this represents the genealogical transmission tree from infectors to infectees.
#' 
#' @importFrom igraph graph
#' @importFrom igraph graph.dfs
#' 
#' @keywords internal
#' 
#' @export
makePhyloTree <- function(Edges, Lengths, Root) {
  G <- graph(edges=t(Edges))
  orderT <- graph.dfs(G,Root)$order
  newLengths=0*Lengths
  
  oE<-order(Edges[,2]) 
  Edges <- Edges[oE, ]  
  Lengths <-Lengths[oE]
  newEdges <- matrix(NA, nrow(Edges), ncol(Edges)) 
  ooT<-order(orderT[-1]);
  newEdges[ooT,] <- Edges
  newLengths[ooT]= Lengths
  # newLengths[newLengths==0]=0.01
  
  Nnode <- (length(orderT)-1)/2;
  Ntips <- Nnode+1; 
  pt <- rtree(Ntips); 
  pt$Nnode <- Nnode; #
  pt$edge <-  newEdges;
  pt$edge.length=newLengths;
  pt$tip.label=paste("t",1:Ntips,sep="")
  return(pt)
}
