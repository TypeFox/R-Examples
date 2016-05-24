makeAnnotList <- function(xAnnot){
  result <- list()
  for(i in 1:nrow(xAnnot)){
    result[[i]] <- xAnnot[i,2:4] 
  }
  names(result) <- xAnnot[,1]
  result
}