## This is the same as the function in quasse2
descendants <- function(node, edge) {
  ans <- node
  repeat {
    node <- edge[edge[,1] %in% node,2]
    if ( length(node) > 0 )
      ans <- c(ans, node)
    else
      break
  }

  unlist(ans)
}

descendants.C <- function(node, edge, n.tip) {
  storage.mode(edge) <- "integer"
  storage.mode(node) <- "integer"
  storage.mode(n.tip) <- "integer"
  .Call("r_descendants", node, edge, n.tip, PACKAGE="diversitree")
}

descendants.flag.C <- function(node, edge, n.tip) {
  storage.mode(edge) <- "integer"
  storage.mode(node) <- "integer"
  storage.mode(n.tip) <- "integer"
  .Call("r_descendants_flag", node, edge, n.tip, PACKAGE="diversitree")
}

descendants.idx.C <- function(node, edge, n.tip) {
  storage.mode(edge) <- "integer"
  storage.mode(node) <- "integer"
  storage.mode(n.tip) <- "integer"
  .Call("r_descendants_idx", node, edge, n.tip, PACKAGE="diversitree")
}

get.children <- function(edge, n.tip) {
  ## To construct the children vector, this is what I used to do:
  ##   lapply(idx[!is.tip], function(x) edge[edge[,1] == x,2])
  ## But this is slow for large trees.  This is faster:
  ## children <- split(edge[,2], edge[,1])
  ## Surprisingly, most of the time is in coercing edge[,1] into a
  ## factor.
  x <- as.integer(edge[,1])
  levels <- as.integer((n.tip+1):max(edge[,1]))
  f <- match(x, levels)
  levels(f) <- as.character(levels)
  class(f) <- "factor"
  children <- split(edge[,2], f)
  names(children) <- NULL

  ## In most cases, this will have been checked by check.tree()
  ## This is currently the time sink here.
  if ( !all(unlist(lapply(children, length)) == 2) )
    stop("Multifircations/unbranched nodes in tree - best get rid of them")

  rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), 2)))
}

ancestors <- function(phy, i=seq_along(phy$tip.label)) {
  anc <- i
  edge <- phy$edge
  while ( any(!is.na(i)) ) {
    i <- edge[match(i, edge[,2]),1]
    anc <- cbind(anc, i, deparse.level=0)
  }

  apply(anc, 1, function(x)
        c(rev(x[!is.na(x)]), rep(NA, sum(is.na(x)))))
}

## Compute the MRCA of tips with indices in 'tips'
mrca.tipset <- function(phy, tips) {
  if ( is.character(tips) )
    tips <- match(tips, phy$tip.label)
  
  if ( length(tips) == 1 )
    tips
  else {
    anc <- ancestors(phy, tips)
    j <- which(apply(anc, 1, function(x) length(unique(x))) > 1)[1]
    anc[j-1,1]
  }
}

## Similar to ape's branching.times(), but returns the height above
## the root node, even for non-ultrametric trees.  Includes tip times.
branching.heights <- function(phy) {
  if (!inherits(phy, "phylo"))
    stop('object "phy" is not of class "phylo"')
  phy <- reorder(phy, "cladewise")

  edge <- phy$edge
  n.node <- phy$Nnode
  n.tip <- length(phy$tip.label)

  ht <- numeric(n.node + n.tip) # zero'd
  for (i in seq_len(nrow(edge)))
    ht[edge[i, 2]] <- ht[edge[i, 1]] + phy$edge.length[i]

  ## Ugly, but fairly compatible with branching.times()
  names.node <- phy$node.label
  if (is.null(names.node))
    names.node <- (n.tip + 1):(n.tip + n.node)
  names(ht) <- c(phy$tip.label, names.node)

  ht
}

## This only works for ultrametric trees:
branching.depth <- function(len, children, order, tips) {
  depth <- numeric(nrow(children))
  depth[tips] <- 0
  for ( i in order )
    depth[i] <- depth[children[i,1]] + len[children[i,1]]
  depth
}

get.descendants <- function(node, tree, tips.only=FALSE,
                            edge.index=FALSE) {
  n.tip <- length(tree$tip.label)
  if ( length(node) != 1 )
    stop("'node' must be scalar")
  if ( is.character(node) ) {
    node <- match(node, tree$node.label) + n.tip
    if ( is.na(node) )
      stop(sprintf("Node '%s' not found in tree"), node)
  } else {
    node <- check.integer(node)
    if ( node >= 1 && node < tree$Nnode ) # on 1..(n.node), probably
      node <- node + n.tip
    else if ( !(node > n.tip && node <= n.tip + tree$Nnode) )
      stop("Invalid node number")
  }
    
  edge <- tree$edge

  desc <- descendants.C(node, edge, n.tip)
  if ( tips.only )
    desc <- desc[desc <= n.tip]
  if ( edge.index )
    desc <- match(desc, edge[,2])
  desc
}
