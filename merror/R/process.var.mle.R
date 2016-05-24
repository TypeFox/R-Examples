"process.var.mle" <-
function (sigma2,s,beta.bars,N,n) 
{
# Jaech p. 186  equations 6.37 - 6.3.10

  d2 <- sum(beta.bars^2/sigma2)
  d3 <- sum(diag(s)/sigma2)
  
  d4.1 <- 0
  
  for(i in 1:N)
  {
    for(j in 1:N)
    {
       if(i!=j) d4.1 <- d4.1 + beta.bars[i]^2*s[j,j]/(sigma2[i]*sigma2[j])
    }
  }
  
  d4.2 <- 0
  
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      d4.2 <- d4.2 + beta.bars[i]*beta.bars[j]*s[i,j]/(sigma2[i]*sigma2[j])
    }
  }
  
  d4 <- d4.1 - 2*d4.2
  
# cat("\nprocess.var.mle sigma.mu2 =",(n-1)*(d2*d3 - d4)/(n*d2^2) - 1/d2)
  
  # return process variance
  (n-1)*(d2*d3 - d4)/(n*d2^2) - 1/d2
}
