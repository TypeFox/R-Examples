UniqueGenotypeCounts <- function(X, verbose = TRUE) {
  n <- nrow(X)
  Ones <- as.matrix(rep(1,nrow(X)),ncol=1)
  X <- cbind(X[,1:3],Ones)
  rownames(X) <- paste("R",1:nrow(X),sep="")
  colnames(X) <- c("AA","AB","BB","O")
  X <- as.data.frame(X)
  Y <- aggregate(X[,4],by=list(X[,1],X[,2],X[,3]),sum)
  colnames(Y) <- c("AA","AB","BB","w")
  nunique <- nrow(Y)
  if(verbose) {
      cat(n,"rows in X\n")
      cat(nunique,"unique rows in X\n")
  }
  return(Y)
}

