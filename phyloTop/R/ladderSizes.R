#' Ladder sizes
#' 
#' Finds the sizes and positions of ladders in the tree. 
#' A ladder is here defined to be a series of consecutive nodes in the tree,
#' each of which has exactly one tip child (as counted by \code{\link{ILnumber}}).
#' The size of the ladder is given by the number of nodes in the chain.
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' 
#' @return A list of:
#' \itemize{
#' \item ladderSizes the sizes of ladders in the tree
#' \item ladderNodes the ladder nodes in the tree
#' \item ladderEdges the edges between ladder nodes of the tree
#' }
#' 
#' @seealso \code{\link{ILnumber}}, \code{\link{ladderShow}}
#'
#' @import ape  
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph cluster.distribution
#' @importFrom igraph components
#' 
#' @examples
#' ## Find ladder sizes in a random tree with 20 tips:
#' tree <- rtree(20)
#' plot(tree)
#' ladderSizes(tree)
#' # note that the ladders can be highlighted in a plot using ladderShow:
#' ladderShow(tree)
#' 
#' 
#' @export
ladderSizes <- function(tree) {
  tree <- phyloCheck(tree)
  ntip=length(tree$tip.label)
  nn=tree$Nnode
  Ancs=(ntip+1):(ntip+nn) # assumes tips are 1:ntip, then internal nodes
  
  # for each internal node (numbered 1:(ntips-1) in the rows of "pointers"), find its immediate children 
  Pointers=t(vapply(Ancs, function(x) tree$edge[tree$edge[,1]==x,2], FUN.VALUE=c(1,2))) 
  
  # create "DP": like pointers, but any tip is replaced with a zero, any internal node with a 1
  DP=Pointers
  DP[DP<=ntip]=0 
  DP[DP>ntip]=1
  
  BrIsLadder=(rowSums(DP)==1) # which branches are ladders
  LadderBr=(1:nrow(DP))[BrIsLadder] # vector of only br that are part of a ladder
  LadderNodes=ntip+LadderBr  # vector of NODES (original node id) that are part of ladders
  allIsLadder=c(rep(0,ntip),BrIsLadder) # 0 for all tips, then 1 or 0 for internal nodes
  
  # make matrix where first column is c(LadderNodes, LadderNodes) and 
  # second column is c(first child of ladder node, second child of ladder node)
  EL=rbind( cbind(LadderNodes, Pointers[LadderBr,1]),cbind(LadderNodes,Pointers[LadderBr,2])) # Ancestor, Descendant where Ancestor is a Ladder Node (but fast as no finding in the edgelist).
  
  ToKeep=allIsLadder[EL[,2]] # vector where entry is 1 only where second col of EL is a ladder node
  EdgeList=EL[ToKeep==1,,drop=FALSE]     # only keep edges where ancestor and des are ladd. drop=FALSE keeps it as a matrix even when it has only one row
  
  if (length(EdgeList)==0) {return(list(ladderSizes=0,ladderNodes=NULL,ladderEdges=NULL))} # no ladder nodes
  
  else {
  chEL <- matrix(as.character(EdgeList),nrow=nrow(EdgeList),ncol=2) # version of EdgeList where entries are characters
  # use igraph to find connected components of this edge list  
  gr=graph_from_edgelist(chEL)
  ladderSizes <- components(gr)$csize # sizes of different ladder components
  
  # list the nodes involved in ladders, in numerical order:
  ladderNodes <- sort(unique(c(EdgeList[,1],EdgeList[,2])))
  
  # find the rows of tree$edge which correspond to rows of edge list - these are the ladder edges
  ladderEdges <- sapply(1:nrow(EdgeList), function(x)
    intersect(which(tree$edge[,1]==EdgeList[x,1]),which(tree$edge[,2]==EdgeList[x,2]))
    )
  
  return(list(ladderSizes=ladderSizes,ladderNodes=ladderNodes,ladderEdges=ladderEdges) )
  }
}