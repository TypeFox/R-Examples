fyshuffle <-
function(x){
  n <- length(x)
  i <- n 
  while (i>0){
    j <- floor(runif(1) * i+1)
    if (i != j){
      temp <- x[i]
      x[i] <- x[j]
      x[j] <- temp
    }
    i <- i - 1
  }
  return(x)
}
