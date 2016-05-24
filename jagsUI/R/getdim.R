#This function gets the dimensions of non-scalar parameters 
#for which the user has requested posterior distributions.

get.dim <- function(params){
  
  #Get all unique parameters (i.e., collapse indexed non-scalars)
  ps <- unique(sapply(strsplit(params, "\\["), "[", 1)) 
  #Slice indexes from non-scalar parameter entries
  test <- sapply(strsplit(params, "\\["), "[", 1)
  
  #Calculate dimension for each parameter
  dim <- lapply(ps, function(i){
    
    w <- params[test==i]
    w2 <- strsplit(w,'\\[')[[length(w)]][2]
    w3 <- strsplit(w2,"\\]")[[1]] 
    as.numeric(unlist(strsplit(w3,",")))
        
  })
  
  names(dim) = ps
  dim
  
}