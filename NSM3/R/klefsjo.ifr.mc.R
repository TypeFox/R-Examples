klefsjo.ifr.mc<-function (x,alternative="two.sided", exact=FALSE, min.reps=100, max.reps=1000, delta=10^-3) 
{
  find.ifr = function(x){
    # IFR
    x.sort = sort(x)
    n = length(x.sort)
    D = n*x.sort[1] #find D_{1}
    n1 = n+1 #use this in calculating alpha
    alpha = ((n1^3) - (3*(n1^2)) + (2*n1))/6
    for(i in 2:n){
      D = c(D,(n-i+1)*(x.sort[i]-x.sort[i-1])) #alpha_{1}
      alpha  = c(alpha, (((n1^3)*i) - (3*(n1^2)*(i^2)) + (2*n1*(i^3)))/6)
    }
    A = sum(alpha*D)/sum(D)
    A.star = A*sqrt(7560/(n^7))
    return(A.star)
  }
  p.g.mc= function (A.star, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.ifr(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn>A.star])/reps
      dsn = c(dsn,find.ifr(rexp(n)))
      if(abs(p-length(dsn[dsn>A.star])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  }
  p.l.mc= function (A.star, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the LESS THAN probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.ifr(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn<A.star])/reps
      dsn = c(dsn,find.ifr(rexp(n)))
      if(abs(p-length(dsn[dsn<A.star])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  } 
  # The next three lines are modified from fisher.test to get the correct alternative
  # hypotheses.  Citation needed?
  alternative <- char.expand(alternative, c("two.sided", "dfr", "ifr"))
  if (length(alternative) > 1L || is.na(alternative)) 
    stop("alternative must be \"two.sided\", \"dfr\" or \"ifr\"")
  
  
  
  A.star = find.ifr(x)
  n = length(x)
  
  # large sample approximation
  if(n>=9){
    if(exact==FALSE){
      
      if(alternative=="dfr"){
        p=pnorm(A.star)
      }
      
      #dmrl
      if(alternative=="ifr"){
        p=pnorm(A.star, lower.tail=F)
      }
      
      #not equal 
      if(alternative=="two.sided"){
        p=2*pnorm(abs(A.star), lower.tail=F)
      }
      
      
      cat("A*=", A.star, "\n", "p=", p,"\n")
      return(list(A.star=A.star,prob=p))
    }
  }
  
  # Exact Test
  #dfr
  if(alternative=="dfr"){
    p=p.l.mc(A.star,n, min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #ifr 
  if(alternative=="ifr"){
    p=p.g.mc(A.star,n,min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(p.l.mc(A.star,n, min.reps=min.reps, max.reps=max.reps, delta=delta),p.g.mc(A.star,n,min.reps=min.reps, max.reps=max.reps, delta=delta))
  }
  
  
  
  cat("A*=", A.star, "\n", "p=", p,"\n")
  return(list(A.star=A.star,p=p))
  
  
}
