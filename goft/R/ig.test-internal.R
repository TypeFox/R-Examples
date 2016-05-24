# Anderson-Darling test for the gamma distribution with shape parameter equal to 0.5
# and unknown scale parameter


`.gammadist.test2` <- function(z)
{      
  n  <- length(z)
  b. <- 2*mean(z)
  s  <- sort(z)/b.  
  theop <-  pgamma(s,.5)
  # Anderson-Darling statistic
  ad <- -n-sum((2*(1:n)-1)*log(theop) + (2*n+1-2*(1:n))*log(1-theop))/n
  m <- .645
  l <- 1.64
  # Approximated p-value
  p.value  <- 1-(pnorm(sqrt(l/ad)*(ad/m-1))+exp(2*l/m)*pnorm(-sqrt(l/ad)*(ad/m+1)))
  results <- list(statistic = c("AD" = ad), p.value =p.value)                  
  class(results) = "htest"
  return(results)  
}

