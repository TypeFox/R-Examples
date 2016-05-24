Binom_Sim <-
function(size,p,N)  {
  q <- 1-p
  x <- numeric(N)
  for(i in 1:N){
    temp <- runif(1)
    j <- 0; cc <- p/(1-p); prob <- (1-p)^size; F <- prob
    while(temp >= F){
      prob <- cc*(size-j)*prob/(j+1); F <- F+prob; j <- j+1
    }
    x[i] <- j
  }
  return(x)
}
