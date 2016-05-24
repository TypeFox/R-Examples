
edgeRetention <- function(matrixAList, threshold, p){
  Ahata <- array( unlist(matrixAList), 
                  dim=c( nrow(matrixAList[[1]]), ncol(matrixAList[[1]]), length(matrixAList)))
  adjacencyMean <- apply(Ahata, c(1,2), function(x) mean(x!=0))
  adjacencyMean[adjacencyMean < threshold] <- 0
  round(100*adjacencyMean)
}
