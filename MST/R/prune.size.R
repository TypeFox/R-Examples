prune.size <-
function(tree){
  if(is.null(dim(tree))) stop("No Need to Prune Further.")
  result <- NULL; n.tmnl <- sum(is.na(tree$var)); subtree <- 1            
  while (n.tmnl > 1 ) {
    # if (n.tmnl==5) {btre <- tree; print(btre)}
    internal <- tree$node[!is.na(tree$cut)]; l <- length(internal); 
    r.value <- 1:l
    for(i in 1:l) {
      branch <- tree[is.element(tree$node,c(internal[i], de(internal[i], tree=tree))),]
      score <- as.numeric(as.vector(branch$score))
      r.value[i] <- sum(score, na.rm=TRUE) / sum(!is.na(score))
    }
    alpha <- min(r.value)
    nod.rm <- internal[r.value == alpha]; 
    # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
    G <- sum(as.numeric(as.vector(tree$score)), na.rm=TRUE);
    G.test <- sum(as.numeric(as.vector(tree$score.test)), na.rm=TRUE)
    result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tree), 
                                  size.tmnl=nrow(tree)-l, alpha=alpha, G=G, G.test=G.test))
    tree <- tree[!is.element(tree$node, de(nod.rm,tree)),]
    tree[match(nod.rm, tree$node), c("var", "vname", "cut", "score", "score.test")] <- NA
    n.tmnl <- sum(is.na(tree$cut))
    subtree <- subtree + 1          
  }
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tree), 
                                size.tmnl=1, alpha=9999, G=0, G.test=0))    
  result <- as.data.frame(result)
  result
}
