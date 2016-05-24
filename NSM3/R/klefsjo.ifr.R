klefsjo.ifr<-function (x,alternative="two.sided", exact=FALSE) 
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
  
  prob.klefsjo = function (t,n) 
  {
    #P(A>t) = sum_{j=1}^{n}(product_{i=1 (i neq j)}^{n}[(a_{j}-t)/(a_{j}-a_{i})]delta_{j}) 
    
    #NOTE: This DOES NOT WORK if two or more a_{j} are equal 
    
    n1 = n+1 #for readability in the formula
    j = 1:n
    coe  = (((n1^3)*j) - (3*(n1^2)*(j^2)) + (2*n1*(j^3)))/6 #calculate the a's
    t = t/(sqrt(7560/(n^7))) # convert from A* to A
    
    
    delta = function(coeff, t){
      #if coeff> t, delta_{j} =1, else 0
      if(coeff>t) return(1)
      else return(0)
    }
    
    sum.val = 0 #initialize sum to 0
    for(j in 1:n){
      prod.val = 1 #initialize product to 1
      for(i in 1:n){
        if(i!=j){
          prod.val = prod.val * (((coe[j]-t)/(coe[j]-coe[i]))*delta(coe[j],t))
          
        }
      }
      sum.val = sum.val + prod.val
    }
    return(sum.val)
    
    
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
    p=1-prob.klefsjo(A.star,n)
  }
  
  #ifr 
  if(alternative=="ifr"){
    p=prob.klefsjo(A.star,n)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(1-prob.klefsjo(A.star,n),prob.klefsjo(A.star,n))
  }
  
  if(is.na(p)){
    print("Large Sample Approximation Used because of Equal Coefficients")
    # Recalculate with large sample approx.  
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
  }
  
  cat("A*=", A.star, "\n", "p=", p,"\n")
  return(list(A.star=A.star,p=p))
  
  
}
