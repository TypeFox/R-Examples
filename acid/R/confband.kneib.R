confband.kneib <-
function(samples, level=0.95)
{
  # at least k functions completely cotained in the band
  n <- nrow(samples)
  p <- ncol(samples)
  
  k <- trunc(level*n)+1
  
  pmean <- apply(samples, 2, mean)
  pqlow <- apply(samples, 2, quantile, probs=(1-level)/2) #quantile for each column, i.e. each point in line
  pqupp <- apply(samples, 2, quantile, probs=1-(1-level)/2)
  
  count <- function(x)
    sum(x>upper) + sum(x<lower)
  
  fac1 <- 1
  lower <- pqlow
  upper <- pqupp
  k1 <- sum(apply(samples, 1, count) == 0)
  
  fac2 <- 1.5
  lower <- pmean + fac2*(pqlow-pmean)
  upper <- pmean + fac2*(pqupp-pmean)
  k2 <- sum(apply(samples, 1, count) == 0)
  
  while(k2 < k)
  {
    fac2 <- 1.5*fac2
    lower <- pmean + fac2*(pqlow-pmean)
    upper <- pmean + fac2*(pqupp-pmean)
    k2 <- sum(apply(samples, 1, count) == 0)
  }
  
  while(k1+1 < k2)
  {
    facnew <- 0.5*(fac1+fac2)
    lower <- pmean + facnew*(pqlow-pmean)
    upper <- pmean + facnew*(pqupp-pmean)
    knew <- sum(apply(samples, 1, count) == 0)
    
    if(knew < k)
    {
      k1 <- knew
      fac1 <- facnew
    }
    if(knew >=k)
    {
      k2 <- knew
      fac2 <- facnew
    }
  }
  
  return(list(lower=lower, upper=upper))
}
