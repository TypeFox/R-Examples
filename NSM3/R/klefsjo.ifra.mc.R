klefsjo.ifra.mc<-function (x, alternative="two.sided", exact=FALSE, min.reps=100, max.reps=1000, delta=10^-3) 
{
  find.ifra = function(x){
    x.sort = sort(x)
    n = length(x.sort)
    D = n*x.sort[1] #find D_{1}
    beta = (-1 + (1-(3*n)-(3*(n^2))) + (2*n) + (3*(n^2)) + (n^3) )/6 #beta_{1}
    for(i in 2:n){
      D = c(D,(n-i+1)*(x.sort[i]-x.sort[i-1])) #alpha_{1}
      beta  = c(beta, ((2*(i^3)) - (3*(i^2)) + i*(1-(3*n)-(3*(n^2))) + (2*n) + (3*(n^2)) + (n^3) )/6)
    }
    B = sum(beta*D)/sum(D)
    B.star = B*sqrt(210/(n^5))
    return(B.star)
  }
  
  p.g.mc= function (B.star, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.ifra(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn>B.star])/reps
      dsn = c(dsn,find.ifra(rexp(n)))
      if(abs(p-length(dsn[dsn>B.star])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  }
  p.l.mc= function (B.star, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the LESS THAN probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.ifra(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn<B.star])/reps
      dsn = c(dsn,find.ifra(rexp(n)))
      if(abs(p-length(dsn[dsn<B.star])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  } 
  # The next three lines are modified from fisher.test to get the correct alternative
  # hypotheses.  Citation needed?
  alternative <- char.expand(alternative, c("two.sided", "dfra", "ifra"))
  if (length(alternative) > 1L || is.na(alternative)) 
    stop("alternative must be \"two.sided\", \"dfra\" or \"ifra\"")
  
  
  B.star = find.ifra(x)
  n = length(x)
  
  # large sample approximation
  if(n>=9){
    if(exact==FALSE){
      
      if(alternative=="dfra"){
        p=pnorm(B.star)
      }
      
      #dmrl
      if(alternative=="ifra"){
        p=pnorm(B.star, lower.tail=F)
      }
      
      #not equal 
      if(alternative=="two.sided"){
        p=2*pnorm(abs(B.star), lower.tail=F)
      }
      
      
      cat("B*=", B.star, "\n", "p=", p,"\n")
      return(list(B.star=B.star,prob=p))
    }
  }
  
  # Exact Test
  #dfr
  if(alternative=="dfra"){
    p=p.l.mc(B.star,n, min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #ifr 
  if(alternative=="ifra"){
    p=p.g.mc(B.star,n,min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(p.l.mc(B.star,n, min.reps=min.reps, max.reps=max.reps, delta=delta),p.g.mc(B.star,n,min.reps=min.reps, max.reps=max.reps, delta=delta))
  }
  
  
  
  cat("B*=", B.star, "\n", "p=", p,"\n")
  return(list(B.star=B.star,p=p))
  
  
}
