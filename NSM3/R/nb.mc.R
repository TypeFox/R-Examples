nb.mc<-function (x, alternative = "two.sided",exact=FALSE, min.reps=100, max.reps=1000, delta=10^-3) 
{
  # 11.1. A Test of Exponentiality versus IFR Alternatives (Epstein)
  # Monte Carlo Version for the exact test.  
  
  # Assumptions: 
  #   A1. x_{i} ~ iid F (F is continuous)
  #   A2. F(a) = 0 for a<0.
  
  # If exact == FALSE, then the large sample approximation will be used if n>=9
  
  find.nbu = function (x) 
    # find's the value of NBU
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
  }
  
  p.g.mc= function (T, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the GREATER THAN probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.nbu(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn>T])/reps
      dsn = c(dsn,find.nbu(rexp(n)))
      if(abs(p-length(dsn[dsn>T])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  }
  
  p.l.mc= function (T, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the LESS THAN probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.nbu(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn<T])/reps
      dsn = c(dsn,find.nbu(rexp(n)))
      if(abs(p-length(dsn[dsn<T])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  }
  
  
  # The next three lines are modified from fisher.test to get the correct alternative
  # hypotheses.  Citation needed?
  alternative <- char.expand(alternative, c("two.sided", "nbu", "nwu"))
  if (length(alternative) > 1L || is.na(alternative)) 
    stop("alternative must be \"two.sided\", \"nbu\" or \"nwu\"")
  
  T = find.nbu(x)
  n = length(x)
  
  # large sample approximation
  if(n>=9){
    if(exact==FALSE){
      b <- T
      e <- n*(n-1)*(n-2)/8
      g <- (3/2)*n*(n-1)*(n-2)
      h <- (5/2592)*(n-3)*(n-4)
      i <- (n-3)*(7/432)
      j <- 1/48
      k <- g*(h+i+j)
      s <- sqrt(k)
      #            p <- 1-pnorm(abs((b-e)/s))
      T.star = (b-e)/s
      #nbu
      if(alternative=="nbu"){
        p=pnorm(T.star)
      }
      
      #nwu
      if(alternative=="nwu"){
        p=pnorm(T.star, lower.tail=F )
      }
      
      #not equal 
      if(alternative=="two.sided"){
        p=2*pnorm(T.star, lower.tail=F)
      }
      cat("T*=", T.star, "\n", "p=", p,"\n")
      return(list(T=T.star,prob=p))
    }
  }
  
  # Exact Test
  #nbu 
  if(alternative=="nbu"){
    p=p.l.mc(T,n, min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #nwu 
  if(alternative=="nwu"){
    p=p.g.mc(T,n,min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(p.l.mc(T,n, min.reps=min.reps, max.reps=max.reps, delta=delta),p.g.mc(T,n,min.reps=min.reps, max.reps=max.reps, delta=delta))
  }
  cat("T=", T, "\n", "p=", p,"\n")
  return(list(T=T,p=p))
  
}
