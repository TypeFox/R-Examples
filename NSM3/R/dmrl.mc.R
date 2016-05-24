dmrl.mc<-function (x,alternative = "two.sided",exact=FALSE, min.reps=100, max.reps=1000, delta=10^-3) 
{
  # A function for 11.3: A test of exponentiality versus
  # DMRL Alternatives (Hollander-Proschan)
  
  find.V.star = function(x){
    x.ord = sort(x)
    n = length(x.ord)
    i = c(1:n)
    c.i = ((4/3)*(i^3)) - (4*n*(i^2)) +(3*(n^2)*i) -
      ((n^3)/2) + ((n^2)/2) - ((i^2)/2) + (i/6)
    V = sum(c.i*x.ord)/(n^4)
    V.star = V/mean(x.ord)
    return(V.star)
  }
  
  p.g.mc= function (V.star, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the GREATER THAN probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.V.star(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn>V.star])/reps
      dsn = c(dsn,find.V.star(rexp(n)))
      if(abs(p-length(dsn[dsn>V.star])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  }
  
  p.l.mc= function (V.star, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the LESS THAN probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.V.star(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn<V.star])/reps
      dsn = c(dsn,find.V.star(rexp(n)))
      if(abs(p-length(dsn[dsn<V.star])/reps)<=delta){
        return(p) #if p converges to be w/i delta, then return
      }
      reps = reps +1
    }
    print("Warning: reached maximum reps without converging within delta")
    return(p)
    
  } 
  # The next three lines are modified from fisher.test to get the correct alternative
  # hypotheses.  Citation needed?
  alternative <- char.expand(alternative, c("two.sided", "dmrl", "imrl"))
  if (length(alternative) > 1L || is.na(alternative)) 
    stop("alternative must be \"two.sided\", \"dmrl\" or \"imrl\"")
  
  V.star = find.V.star(x)
  n = length(x)
  
  # large sample approximation
  if(n>=9){
    if(exact==FALSE){
      V.prime = sqrt(210*n)*V.star
      
      #imrl
      if(alternative=="imrl"){
        p=pnorm(V.prime)
      }
      
      #dmrl
      if(alternative=="dmrl"){
        p=pnorm(V.prime, lower.tail=F )
      }
      
      #not equal 
      if(alternative=="two.sided"){
        p=2*pnorm(V.prime, lower.tail=F)
      }
      cat("V'=", V.prime, "\n", "p=", p,"\n")
      return(list(V=V.prime,prob=p))
    }
  }
  
  # Exact Test
  #imrl 
  if(alternative=="imrl"){
    p=p.l.mc(V.star,n, min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #dmrl 
  if(alternative=="dmrl"){
    p=p.g.mc(V.star,n,min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(p.l.mc(V.star,n, min.reps=min.reps, max.reps=max.reps, delta=delta),p.g.mc(V.star,n,min.reps=min.reps, max.reps=max.reps, delta=delta))
  }
  cat("V*=", V.star, "\n", "p=", p,"\n")
  return(list(V=V.star,p=p))
  
  
  
  
}
