"quantregForest" <-
function(rf, ...){

  ## Some checks 
  if(! class(rf) %in% c("randomForest") )
    stop(" rf must be of class `randomForest' ")
  
  
  nodesX <- attr(predict(rf,x,nodes=TRUE),"nodes")
  rownames(nodesX) <- NULL
  nnodes <- max(nodesX)
  ntree <- ncol(nodesX)
  n <- nrow(x)
  valuesNodes <- matrix(nrow=nnodes,ncol=ntree)
  for (tree in 1:ntree){
      shuffledNodes <- nodesX[rank(ind <- sample(1:n,n)),tree]
      useNodes <- sort(unique(as.numeric(shuffledNodes)))
      valuesNodes[useNodes,tree] <- y[ind[match(useNodes,shuffledNodes )]]
  }
  
  class(rf) <- c("quantregForest","randomForest")

  qrf[["valuesNodes"]] <- valuesNodes

  
  return(qrf)
}
