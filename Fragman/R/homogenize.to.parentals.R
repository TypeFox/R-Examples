homogenize.to.parentals <-
function(x, parents, window){
  #d1 <- data.frame(matrix(parents, nrow=1))
  
  alls <- numeric()
  for(i in 1:length(x$wei)){
    bb <- abs(x$wei[i] - parents)
    bb2 <- which(bb < window)[1]
    if(is.na(bb2)){
      alls[i] <- x$wei[1]
    }else{alls[i] <- parents[[bb2]]}
    
  }
  x$wei <- alls
  # returns alleles but in parental conformation
  return(x)
}
