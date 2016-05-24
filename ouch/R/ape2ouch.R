# This function takes a tree file in the phylo format (from
# ape/read.tree) and creates an ouchtree object.
# The option exists whereby users can change branch lengths
# while keeping the same topology.
#
#  t = object of type 'phylo', as returned by ape/read.tree
#  branch.lengths = optional branch length vector in same
#   order as t$edge.length; default is to use t$edge.length
# 
# D. Ackerly, dackerly at berkeley dot edu, July 18, 2006.
# Modified by A.A. King, kingaa at umich dot edu, 5/15/2007.

ape2ouch <- function (tree, scale = TRUE, branch.lengths = tree$edge.length) {

  if (!inherits(tree,'phylo'))
    stop(sQuote("tree")," must be of class ",sQuote("phylo"))
  
  nnodes <- nrow(tree$edge)+1              # number of nodes
  n.term <- length(tree$tip.label)         # number of terminal nodes
  n.int <- nnodes-n.term                # internal nodes

  tmp <- matrix(NA,nnodes,2)
  tmp[-1,1:2] <- tree$edge
  tmp[1,2] <- nnodes+1
  bl <- c(0,branch.lengths)
  reord <- order(-tmp[,2])
  bl <- bl[reord]
  tmp <- tmp[reord,]

  node <- seq(nnodes)
  ancestor <- rep(NA,nnodes)

  for (n in 2:nnodes) {
    anc <- which(tmp[,2]==tmp[n,1])
    if (length(anc)>1) stop("invalid tree")
    if (length(anc)>0) {            # the node has a non-root ancestor
      ancestor[n] <- node[anc]
    } else {                    # the node has the root as an ancestor
      ancestor[n] <- node[1]
    }
  }
  
  if (is.null(tree$node.label))
    tree$node.label <- rep('',n.int)
  species <- rev(c(tree$tip.label,tree$node.label[-1],tree$node.label[1]))

  times <- rep(NA,nnodes)
  for (n in 1:nnodes)
    times[n] <- branch.height(node,ancestor,bl,n)
  if (is.logical(scale)) {
    if (is.na(scale)) stop("if ",sQuote("scale")," is logical, it must be either true or false")
    if (scale) times <- times/max(times)
  } else if (is.numeric(scale)) {
    times <- times/abs(scale)
  } else {
    stop(sQuote("scale")," must be either logical or numeric")
  }
  
  ouchtree(
           nodes=node,
           ancestors=ancestor,
           times=times,
           labels=species
           )
}

branch.height <- function (node, anc, bl, k) {
  ## this recursion might prove expensive for large trees
  ## there should be a direct method 
  if (is.na(anc[k])) {
    0
  } else {
    bl[k]+Recall(node,anc,bl,anc[k])
  }
}
