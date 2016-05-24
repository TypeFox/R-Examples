obtain.btree <-
function(tree, bsize=6){
  btre <- NULL
  if (is.factor(bsize)) bsize <- as.numeric.factor(bsize)
  if (bsize==1) { btre <- tree[1,]; btre[, c("var", "cut", "score", "score.test")] <- NA}
  else if (bsize <1) stop("THE BEST TREE SIZE bsize= MUST BE >=1!")
  else {
    n.tmnl <- sum(is.na(tree$cut)); indicator <- TRUE
    if (bsize > n.tmnl){stop("THE BEST TREE SIZE bsize PROVIDED IS LARGER THAN THE FULL TREE THAT YOU HAVE GROWN.")
    } else if (bsize == n.tmnl){btre <- tree;  indicator <- FALSE}
    while (n.tmnl >= bsize && indicator ==TRUE) {
      # print(tree); print(cbind(n.tmnl, bsize))
      internal <- tree$node[!is.na(tree$cut)]; l <- length(internal); 
      r.value <- 1:l
      for(i in 1:l) {
        branch <- tree[is.element(tree$node,c(internal[i], de(internal[i], tree=tree))),]
        score <- as.numeric(as.vector(branch$score))
        r.value[i] <- sum(score, na.rm=TRUE) / sum(!is.na(score))
      }
      alpha <- min(r.value)
      nod.rm <- internal[r.value == alpha]; 
      tree <- tree[!is.element(tree$node, de(nod.rm,tree)),]
      tree[match(nod.rm, tree$node), c("var", "vname", "cut", "score", "score.test")] <- NA
      n.tmnl <- sum(is.na(tree$cut))
      # print(cbind(n.tmnl, bsize)) 
      if (n.tmnl==bsize) {btre <- tree;  indicator <- FALSE}
    }
  }
  if (is.null(btre)) print(paste("The optimally-pruned subtree sequence does not have a subtree of bsize = ", bsize, sep="")) 
  return(btre)
}
