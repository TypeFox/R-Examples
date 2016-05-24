"sigma_mle" <-
function (i,s,sigma2,sigma.mu2,beta.bars,N,n) 
{
# Jaech p. 185-186 equations 6.3.1 - 6.3.6

# this is the ith iteration

# if(sigma.mu2 < 0) sigma.mu2 <- 0

  a <- beta.bars*sqrt(sigma.mu2)
  
#  cat("\na=",a)
  
  b0 <- a[i]^2
  
#  cat("\nb0=",b0)
  
  b1 <- 0
  for(j in 1:N) { if(j!=i) { b1 <- b1 + a[j]^2/sigma2[j] }}
  b1 <- b1 + 1
  
#  cat("\nb1=",b1)
  
  b2.1 <- 0
  for(j in 1:N) { if(j!=i) { b2.1 <- b2.1 + s[j,j]/sigma2[j] }}
  
#  cat("\nb2.1=",b2.1)
  
  b2.2 <- 0
  for(j in 1:N) { if(j!=i) { b2.2 <- b2.2 + a[j]*s[i,j]/sigma2[j] }}
  
#  cat("\nb2.2=",b2.2)
  
  b2 <- b1*s[i,i] + a[i]^2*b2.1 - 2*a[i]*b2.2
  
#  cat("\nb2=",b2)
  
  b3.1 <-0
  for(j in 1:N) { if(j!=i) { b3.1 <- b3.1 + s[j,j]/sigma2[j] }}
  
#  cat("\nb3.1=",b3.1)

  b3.2 <-0
  for(j in 1:N) { if(j!=i) { b3.2 <- b3.2 + a[j]^2*s[j,j]/sigma2[j]^2}}
  
#  cat("\nb3.2=",b3.2)
  
  b3.3 <-0
  for(j in 1:N)
  {
    if(j!=i) 
    {
      for(k in 1:N)
      {
        if(k>j & k!=i) b3.3 <- b3.3 + a[j]*a[k]*s[j,k]/(sigma2[j]*sigma2[k])
      }
    }  
  }
  b3 <- b1*b3.1 - b3.2 - 2*b3.3
  
#  cat("\nb3=",b3)
  
  (n-1)*(b1*b2 - b0*b3)/(n*b1^2) - b0/b1 # return sigma2[i]     

}
