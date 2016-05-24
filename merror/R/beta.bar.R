"beta.bar" <-
function (x) 
{
# compute beta.bar (Jaech, p. 184)
# model for ith instrument and kth item:
# x[i,k] = alpha[i] + beta[i]mu[k] + epsilon[i,k]

  # no. of instruments
  N <- dim(x)[2]
  denominator <- vector("numeric",N)
  
  s <- var(x)
 
  numerator <- (apply(s,1,prod)/diag(s))^(1/N)
  
  for(i in 1:N)
  {
    denominator[i] <- prod(s[-i,-i][upper.tri(s[-i,-i])])^(2/(N*(N-2)))
  }
  
  numerator/denominator # beta[i]
  
}
