rFrechet <- function(n){
  -1/log(runif(n))
}

rMaxAR <- function(n,theta)
{
  x <- rFrechet(n)
  
  y <- rbind(rep(x[1],2),cbind((1-theta)*x[-n], x[-1]))
  
  y <- apply(y,1, max)
  
  y

}
