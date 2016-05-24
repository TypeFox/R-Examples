target.help <-
function(genes) {
  T <- diag(1,length(genes))                              
  for (i in 2:length(genes)) {                            
   for(j in 1:(i-1)) {                                    
    T[j,i] <- T[i,j] <- check.path(genes[[i]],genes[[j]]) 
   }
  }
  return(T)
}

