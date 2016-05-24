

#This function takes boxes of the format output by sdprim and converts them
#to the format from original box testing, to take advantage of other
#pre-existing functions


boxconverter <- function(sdbox){

#  cat("WARNING! This function currently only tracks dimension index and not","\n")
#  cat("variable names, which means it will mess up if the box was generated on data", "\n")
#  cat("with a different number of or different ordering of inputs.","\n")

#  sdbox$dimlist$lower
  
#  sdbox$dimlist$upper
  
#  sdbox$box
  
  dimvect <- c(1:length(sdbox$dimlist$either))[sdbox$dimlist$either]  #naive - assumes columns all there and in right order
  
  dimmat  <- matrix(Inf, nrow=sum(sdbox$dimlist$either),ncol=2)
  
  dimmat[,1] <- -Inf
  
  for (i in 1:length(dimvect)){
  
    if (sdbox$dimlist$lower[dimvect[i]]){
  
      dimmat[i,1] <- sdbox$box[1,dimvect[i]]
      
    } 
    
    if (sdbox$dimlist$upper[dimvect[i]]){
  
      dimmat[i,2] <- sdbox$box[2,dimvect[i]]
      
    }
    
  }
  
  newbox <- list(dimvect,dimmat) 
  
  return(newbox)
  
}

  
    
    
  
  