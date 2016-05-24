
edgeSelection <- function(matrixA, q,  nodewise=FALSE){
  # initialize matrix
  p <- ncol(matrixA)
  AhatSelected <- matrix(0, p, p)

  quse <- sample( floor(c(q,q+1)), 1,prob= c(1-q+floor(q),q-floor(q)))
      
  if(!nodewise){
      
      ## choose numberOfEdges largest ones
      selected <- sort(abs(matrixA), decreasing =  TRUE)[1:quse]
      indicesSelected <- which(abs(matrixA) >= min(selected), arr.ind = TRUE)
      
      ## choose elements
      AhatSelected[indicesSelected] <- 1
  }else{
      for (i in 1:nrow(matrixA)){
        AhatSelected[ i,order(abs(matrixA)[i,],rnorm(p),decreasing=TRUE)[1:quse]] <- 1
        AhatSelected[ order(abs(matrixA)[,i],rnorm(p),decreasing=TRUE)[1:quse],i] <- 1
      }
  }
  
  ## return matrix
  AhatSelected
}
