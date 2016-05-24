epstein<-function (x, alternative = "two.sided",exact=FALSE) 
{
  # 11.1. A Test of Exponentiality versus IFR Alternatives (Epstein)
  
  # Assumptions: 
  #   A1. x_{i} ~ iid F (F is continuous)
  #   A2. F(a) = 0 for a<0.
  
  # If exact == FALSE, then the large sample approximation will be used if n>=9
  
  
  p.sum.unif= function (x,n) 
  {
    # Finds the probability of the sum of n uniform(0,1) r.v.'s
    # Formula from An Introduction to Probability Theory and It's Applications by William Feller volume 1 ed. 3 (pg 285)
    # (1/n!)Sum_{i=0}^{n}[(-1^i)(nCi)(x-i)^n] with 0<=i<x
    s = 0
    for(i in 0:n){
      s = s + ((-1)^i)*choose(n,i)*max(0,(x-i))^n
    }
    return(s/factorial(n))
    
  }
  
  
  # The next three lines are modified from fisher.test to get the correct alternative
  # hypotheses.  Citation needed?
  alternative <- char.expand(alternative, c("two.sided", "dfr", "ifr"))
  if (length(alternative) > 1L || is.na(alternative)) 
    stop("alternative must be \"two.sided\", \"dfr\" or \"ifr\"")
  
  
  x.ord = sort(x) #ordered values
  n = length(x)
  
  D = n*x.ord[1] # D_{1}
  for(i in 2:n){
    D = c(D, (n-i+1)*(x.ord[i]-x.ord[i-1]))
  }
  S = cumsum(D)
  E = sum(S[1:(n-1)])/S[n] #E is script E in the text
  
  
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
    p=1-p.sum.unif(((n-1)/2)-E,n-1)
  }
  
  #ifr 
  if(alternative=="ifr"){
    p=1-p.sum.unif(E,n-1)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(1-p.sum.unif(((n-1)/2)-E,n-1),1-p.sum.unif(E,n-1))
  }
  cat("E=", E, "\n", "p=", p,"\n")
  return(list(E=E,p=p))
  
}
