"precision.grubbs.cb.pd" <-
function(x)
{
# Grubbs ME precision estimates for Constant Bias Model using Paired Data
# Confidence intervals not computed (because we should use MLE instead!)
#
# Jaech 1985, Chapters 3 & 4, p. 144 in particular
#
# Written for any no. of methods N
#
# x[i,k] = alpha[i] + beta[i]*mu[k] + epsilon[i,k]
#   with beta[1] = beta[2] = ... = beta[N]
#   N = no. of methods
#   n = no. of items
# 
# ME (method of moments estimator) and MLE are the same for i=3 instruments
#  except for a factor of (n-1)/n: MLE = (n-1)/n * ME
#  so to get the MLE just multiply the output by (n-1)/n

  n <- dim(x)[1]
  N <- dim(x)[2]

# compute paired differences 

  y <- array(NA,c(N,N,n)) # not all elements are used
  
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      y[i,j,] <- x[,i] - x[,j]
    }
  }
  
# print(y[1,2,])

# compute variances for each paired difference

  V <- matrix(NA,N,N) # not all elements are used
  
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      V[i,j] <- var(y[i,j,])
      V[j,i] <- var(y[i,j,])
    }
  }
  
#  print(V)

# Compute Grubbs (Squared) Precions estimates - stored in r

  S <- vector("numeric",N)
  
  for(i in 1:N)
  {
    S[i] <- 0
    for(j in 1:N)
    {
      if(j!=i) S[i] <- S[i] + V[i,j]
    }
  }
  
# print(S)
  
  VT <- sum(S)/2  

  r <- vector("numeric",N)
  
  for(i in 1:N)
  {
    r[i] <- ((N-1)*S[i] - VT)/((N-1)*(N-2))
  }

#  print(r)

# Return squared precisons
  r 
}
