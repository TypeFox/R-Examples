leafset.prob2 <- function(leafset, y) {
  p <- y$p
  tree <- y$tree
  leafset <- unique(leafset)
  pos <- match(leafset, y$var.names)
  parent <- match(tree[1, ], tree[2, ])
  ##indp, indq correspond to the sets E_{j1} and E_{j2} in the paper
  indp <- indq <- rep(0, ncol(tree))
  activenodes <- pos
  while (any(activenodes %in% tree[2, ])) {
    edges <- match(activenodes, tree[2, ])
    edges <- edges[!is.na(edges)]
    indp[edges] <- 1
    activenodes <- unique(tree[1, edges])
  }
  for (k in 1:ncol(tree)) {
    if (indp[k] == 0 & (is.na(parent[k]) | indp[parent[k]] == 1)) 
      indq[k] <- 1
  }
  ##compute q out of p:
  q <- 1 - p
  succ <- lapply(1:ncol(tree), function(x) which(tree[1, ] == tree[2, x]))
  waiting <- sapply(succ, length) > 0
  while(any(waiting)) {
    Next <- which(waiting & sapply(succ, function(x)(!any(waiting[x]))))
    q[Next] <- q[Next] + p[Next]*sapply(succ[Next], function(x)prod(q[x]))
    waiting[Next] <- FALSE
  }
  return(prod(p[indp==1])*prod(q[indq==1]))
}

