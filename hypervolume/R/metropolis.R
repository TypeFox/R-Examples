#Define the sampling function
rwmetro <- function(target,N,x,burnin,delta)
# target - function to estimate
# N - number of samples requested
# x - starting point for the algorithm
{
  VCOV = delta * diag(length(x))
  
  samples <- matrix(NA, nrow=N+burnin,ncol=length(x))
  for (i in 1:(burnin+N))
  {   
    prop <- MASS::mvrnorm(n = 1, x, VCOV)
    if (runif(1) < min(1, target(prop)/target(x))) 
    {
      x <- prop
    }
    samples[i,] <- x
    
    if (i %% 1000 == 0)
    {
      print(i-burnin)
    }
  }
  return(samples[(burnin):(N+burnin),])
}