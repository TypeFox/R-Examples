#' Find the depth of each node
#' 
#' Determines the depth of each node, defined as the number of steps from the root.
#' (So the root has depth zero, its children have depth one, etc.)
#' The output is given as a list of two vectors: \code{tipDepths} gives the depths of the tips, and \code{nodeDepths} gives the depths of the internal nodes.
#' This replaces the deprecated \code{dists} function.
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param tree a tree of class \code{phylo} or \code{phylo4}. The tree should be binary and rooted; if not it will be coerced into a binary rooted tree using multi2di, if possible.
#' @return A list of two vectors: \code{tipDepths} gives the depths of the tips, and \code{nodeDepths} gives the depths of the internal nodes.
#' 
#' @seealso \code{\link{nodeDepth}}, \code{\link{nodeDepthFrac}}
#' 
#' @import ape
#'   
#' @examples
#' ## Find the node depths in a random tree:
#' tree <- rtree(20)
#' treeDepths <- getDepths(tree)
#' ## Plot tree with node depths labelled:
#' tree$tip.label <- treeDepths$tipDepths
#' tree$node.label <- treeDepths$nodeDepths
#' plot(tree, show.node.label=TRUE)
#' 
#' 
#' @export
getDepths=function(tree) {
  # perform tree checks:
  tree <- phyloCheck(tree)
  
  # reorder (if necessary) edges into helpful order, with the root first
  # Parents appear before descendants
  edge_order <- ape::reorder.phylo(tree, "cladewise", index.only=T)
  edges <- tree$edge[edge_order,]
  edge_lengths <- tree$edge.length[edge_order]
  root <- edges[1,1]
  
  ntip <- length(tree$tip.label)
  nn <- tree$Nnode
  
  # list of ancestral (internal) node numbers
  Ancs <- (ntip+1):(ntip+nn)
  
  # Pointers is a matrix where rows correspond to internal nodes in the order (ntip+1, ..., ntip+nn)
  # The two columns give the two descendants of that node
  Pointers=t(vapply(Ancs, function(x) tree$edge[tree$edge[,1]==x,2], FUN.VALUE=c(1,2))) 
  
  # initialise pDep, the node depth array
  pDep <- NA + Pointers
  # fill in first value: the root is depth 1
  pDep[Ancs==root,] <- c(1,1)
  
  # get the root descendants
  toDo <- Pointers[Ancs==root,]
  # retain any which are not tips - these will be the next entries "to do" in pDep:
  toDo <- toDo[toDo>ntip]
  sdep <- 1
  while (any(is.na(pDep))) {
    # fill in the "to do" entries of pDep with the current depth
    pDep[toDo-ntip,] <- sdep + 1
    # establish their non-tip descendants, as the next "to do" in pDep
    toDo <- Pointers[toDo-ntip,][Pointers[toDo-ntip,]>ntip]
    sdep <- sdep + 1
  }
  oP <- order(Pointers)
  tipDepths <- pDep[oP][1:ntip]
  nodeDepths <- c(0,pDep[oP][Ancs[-length(Ancs)]])
  return(list(tipDepths=tipDepths, nodeDepths=nodeDepths))
}