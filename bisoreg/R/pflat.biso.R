`pflat.biso` <- function(obj,...){
  M <- obj$m
  draws <- obj$postdraws
  mean(rowSums(draws[,1:M])==0)
  }

