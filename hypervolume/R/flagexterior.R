expectednumber <- function(density, rmin, n)
{
  pi^(n/2) * rmin^n / gamma(n/2+1) * density
}

flagexterior <- function(input, density, fuzzfactor)
{  
  distances <- as.matrix(dist(input))
  
  n <- ncol(input)
  
  criticaldistance <- fuzzfactor * (density)^(-1/n)
  
  finalstat <- rep(NA, nrow(input))
  
  for (i in 1:nrow(input))
  {
    # find all the neighbors in the critical ball
    nnids <- which(distances[i,] < criticaldistance & distances[i,] > 0)
    
    # compare to expected number of points
    finalstat[i] <- length(nnids) / expectednumber(density, rmin=criticaldistance, n=n)
    
    print(i)
  }

  
  return(finalstat)
}

#td <- matrix(runif(10000),nrow=500,ncol=2)

#distances <- as.matrix(dist(td))


#fe <- flagexterior(td, density=10000, fuzzfactor=2); plot(td, col=ifelse(fe < 0.1,'red','blue'))