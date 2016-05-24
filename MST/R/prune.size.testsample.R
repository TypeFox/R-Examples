prune.size.testsample <-
function(tree){
  out <- as.list(NULL)
  ntest <- as.numeric(tree[1, ncol(tree)])
  if(is.null(dim(tree))) stop("No Need to Prune Further.")
  result <- NULL; n.tmnl <- sum(is.na(tree$var)); subtree <- 1
  a <- cbind(Ga.2=2, Ga.3=3, Ga.4=4, Ga.log_n=log(ntest))
  max.Ga <- rep(-1e20, 4); size <- rep(0, 4); btree <-as.list(1:4) 
  while (n.tmnl > 1 ) {
    # print(tree)
    internal <- tree$node[!is.na(tree$cut)]; l <- length(internal); 
    r.value <- 1:l
    for(i in 1:l) {
      branch <- tree[is.element(tree$node,c(internal[i], de(internal[i], tree=tree))),]
      score <- as.numeric(as.vector(branch$score))
      r.value[i] <- sum(score, na.rm=TRUE) / sum(!is.na(score))
    }
    alpha <- min(r.value)
    nod.rm <- internal[r.value == alpha];
    G <- sum(as.numeric(as.vector(tree$score.test)), na.rm=TRUE); 
    Ga <- G - a*l
    for (k in 1:4){if (Ga[k] > max.Ga[k]) {max.Ga[k] <- Ga[k]; size[k] <- n.tmnl; btree[[k]] <- tree}}
    result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tree),
                                  size.tmnl=nrow(tree)-l, alpha=alpha, G=G, Ga))
    tree <- tree[!is.element(tree$node, de(nod.rm,tree)),]
    tree[match(nod.rm, tree$node), c("var", "vname", "cut", "score", "score.test")] <- NA
    n.tmnl <- sum(is.na(tree$cut))
    if (n.tmnl ==1) {
      for (k in 1:4){
        if (0 > max.Ga[k]) {max.Ga[k] <- 0; size[k] <- 1; btree[[k]] <- tree}
      }
    }
    subtree <- subtree + 1
  }
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tree),
                                size.tmnl=1, alpha=9999, G=0, Ga=cbind(Ga.2=0, Ga.3=0, Ga.4=0, Ga.log_n=0)))
  result <- as.data.frame(result)
  out$result <- result; out$size <- size; out$btree <- btree
  out 
}
