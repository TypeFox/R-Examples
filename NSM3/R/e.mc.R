e.mc<-function (x, alternative = "two.sided",exact=FALSE, min.reps=100, max.reps=1000, delta=10^-3) 
{
  # 11.1. A Test of Exponentiality versus IFR Alternatives (Epstein)
  # Monte Carlo Version for the exact test.  
  
  # Assumptions: 
  #   A1. x_{i} ~ iid F (F is continuous)
  #   A2. F(a) = 0 for a<0.
  
  # If exact == FALSE, then the large sample approximation will be used if n>=9
  
  find.Epstein = function (x) 
    # find's the value of Epstein's E
  {
    x.ord = sort(x) #ordered values
    n = length(x)
    
    D = n*x.ord[1] # D_{1}
    for(i in 2:n){
      D = c(D, (n-i+1)*(x.ord[i]-x.ord[i-1]))
    }
    S = cumsum(D)
    E = sum(S[1:(n-1)])/S[n] #E is script E in the text
    return(E)
  }
  
  p.mc= function (E, n, min.reps=100, max.reps=1000, delta=10^-3) 
  {
    # returns the monte carlo estimate for the probability.  
    dsn = numeric() #initialize
    for(i in 1:min.reps){
      dsn = c(dsn,find.Epstein(rexp(n)))
    }
    reps = min.reps
    while(reps<=max.reps){
      p = length(dsn[dsn>E])/reps
      dsn = c(dsn,find.Epstein(rexp(n)))
      if(abs(p-length(dsn[dsn>E])/reps)<=delta){
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
  
  E = find.Epstein(x)
  n = length(x)
  
  # large sample approximation
  if(n>=9){
    if(exact==FALSE){
      E.star = (E-((n-1)/2))/sqrt((n-1)/12)
      
      #dfr
      if(alternative=="dfr"){
        p=pnorm(E.star)
      }
      
      #ifr
      if(alternative=="ifr"){
        p=pnorm(E.star, lower.tail=F)
      }
      
      #not equal 
      if(alternative=="two.sided"){
        p=2*pnorm(E.star, lower.tail=F)
      }
      cat("E*=", E.star, "\n", "p=", p,"\n")
      return(list(E=E.star,prob=p))
    }
  }
  
  # Exact Test
  #dfr 
  if(alternative=="dfr"){
    p=p.mc(((n-1)/2)-E,n-1, min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #ifr 
  if(alternative=="ifr"){
    p=p.mc(E,n-1,min.reps=min.reps, max.reps=max.reps, delta=delta)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(p.mc(((n-1)/2)-E,n-1,min.reps=min.reps, max.reps=max.reps, delta=delta),
            p.mc(E,n-1,min.reps=min.reps, max.reps=max.reps, delta=delta))
  }
  cat("E=", E, "\n", "p=", p,"\n")
  return(list(E=E,p=p))
  
}
