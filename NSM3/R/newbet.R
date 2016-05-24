newbet<-function (x) 
{
  n=length(x)
  y = sort(x)
  a = array(0,c(n,n,n))
  
  #computes the NBU statistic
  for(i in 3:n){
    for(j in 2:(i-1)){
      for(k in 1:(j-1)){
        if(y[i]>y[j]+y[k])
          a[i,j,k]<-1
        else if (y[i]==y[j]+y[k])
          a[i,j,k]<-(1/2)
        else a[i,j,k]<-0
      }
    }
  }
  sum(a)
  
  #computes the one-sided P-value
  b <- sum(a)
  e <- n*(n-1)*(n-2)/8
  g <- (3/2)*n*(n-1)*(n-2)
  h <- (5/2592)*(n-3)*(n-4)
  i <- (n-3)*(7/432)
  j <- 1/48
  k <- g*(h+i+j)
  s <- sqrt(k)
  p <- 1-pnorm(abs((b-e)/s))
  results <- list(sum(a), (b-e)/s,p)
  results
}
