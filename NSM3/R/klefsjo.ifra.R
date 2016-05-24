klefsjo.ifra<-function (x, alternative="two.sided", exact=FALSE) 
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
  
  prob.klefsjo=function (t,n) 
  {
    #P(A>t) = sum_{j=1}^{n}(product_{i=1 (i neq j)}^{n}[(a_{j}-t)/(a_{j}-a_{i})]delta_{j}) 
    
    #NOTE:  This works for P(B>t), just change all the a's (alpha in the text) to b's (beta in the text)
    
    #NOTE: This DOES NOT WORK if two or more b_{j} are equal 
    
    j = 1:n
    #calculate the b's
    coe  = ((2*(j^3)) - (3*(j^2)) + j*(1-(3*n)-(3*(n^2))) + (2*n) + (3*(n^2)) + (n^3) )/6
    t = t/(sqrt(210/(n^5))) #convert from B* to B 
    
    
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
    p= 1-prob.klefsjo(B.star,n)
  }
  
  #ifr 
  if(alternative=="ifra"){
    p=prob.klefsjo(B.star,n)
  }
  
  #not equal 
  if(alternative=="two.sided"){
    p=2*min(1-prob.klefsjo(B.star,n),prob.klefsjo(B.star,n))
  }
  
  if(is.na(p)){
    print("Large Sample Approximation Used because of Equal Coefficients")
    # Recalculate with large sample approx.  
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
  }
  
  
  cat("B*=", B.star, "\n", "p=", p,"\n")
  return(list(B.star=B.star,p=p))
  
  
}
