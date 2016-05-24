combpaste <-
function(x,y)
# combpaste()  -  paste elements of x anyd y together in all combinations
{
  n <- 0
  combxy <- rep(0,length(x) * length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      n <- n+1
      combxy[n] <- paste(x[i],":",y[j],sep="",collapse=NULL)
    }
  }
  return(combxy)
}
