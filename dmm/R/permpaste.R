permpaste <-
function(x)
# permpaste() - paste elements of vector x in all permutations pairwise
{
  n <- 0
  permx <- rep(0,length(x)^2)
  for(i in 1: length(x)){
    for(j in 1: length(x)){
      n <- n + 1
      permx[n] <- paste(x[i],":",x[j],sep="", collapse=NULL)
    }
  }
  return(permx)
}
