as.stratifiedvegdata<-function(X) {
  if(!inherits(X,"list")) stop("Input should be of class 'list'")
  n = length(X)
  if(n==0) stop("The list should be of positive length")
  for(i in 1:n) if(!inherits(X[[i]],"matrix") && !inherits(X[[i]],"data.frame")) stop("The list should contain 'data.frame' or 'matrix' objects only")
  p = nrow(X[[1]])
  s = ncol(X[[1]])
  for(i in 2:n) {
    if(nrow(X[[i]])!=p) stop("All elements should have the same number of rows (species)")
    if(ncol(X[[i]])!=s) stop("All elements should have the same number of columns (strata)")
  }  
  class(X)<-c("list","stratifiedvegdata")
  if(is.null(names(X))) names(X)<-1:length(X)
  return(X)
}