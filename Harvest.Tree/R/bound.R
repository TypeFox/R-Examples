
bound <- function(noden)
{
  boundr <- rep(NA,nrow(noden$bounds))
  for(i in 1:nrow(noden$bounds)){
    if(!is.na(noden$bounds[i,3])){
      boundr[i] <- 1
    } 
    else 
      boundr[i] <- sum(!is.infinite(noden$bounds[i,1]),!is.infinite(noden$bounds[i,2]))
  }
  bound <- boundr[boundr!=0]
  return(bound)
}


