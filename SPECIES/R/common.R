##Poisson mixture density
untrunPmix= function(i,p,pi){
  ## =============================================================================
  ## Purpose: This function calculates Poisson mixture density
  ## input:   i  --- integer or integer vector (>=0)
  ##          p  --- mean parameters for Poisson mixture components
  ##          pi --- weight parameters for Poisson mixture components
  ##            
  ## output: the Poisson mixture porbabilities
  ## =============================================================================

  prob = numeric(length(i))
  for (k in 1:length(i)){
    prob[k]=sum(pi[p>0]*dpois(i[k],p[p>0]))
  }
  prob
}

##Zero-truncated Poisson  mixture density
Pmix= function(i,p,pi){

  ## ===============================================================================
  ## Purpose: This function calculates zero-truncated Poisson mixture density
  ## input:   i  --- integer or integer vector (>0)
  ##          p  --- mean parameters for zero-truncated Poisson mixture components
  ##          pi --- weight parameters for zero-truncated Poisson mixture components
  ##            
  ## output: the zero-truncated Poisson mixture porbabilities
  ## =============================================================================

  prob = numeric(length(i))
  for (k in 1:length(i)){
    prob[k]=sum(pi[p>0]*dpois(i[k],p[p>0])/(1-exp(-p[p>0])))
  }	      
  prob
}

##Poisson Gamma mixture density calculation funciton

PmixSCon= function(i,p,pi,a){

  ## ===============================================================================
  ## Purpose: This function calculates zero-truncated Poisson compound Gamma density
  ## input:   i  --- integer or integer vector (>0)
  ##          p  --- mean parameters for Gamma components in the compound Gamma model
  ##          pi --- weight parameters for Gamma components
  ##            
  ## output: the zero-truncated Poisson compound Gamma porbabilities
  ## =============================================================================

  prob = numeric(length(i))
  for (k in 1:length(i)){
    prob[k]=sum(pi[p>0]*exp((a+i[k]-2)*log(2)+log(dbeta(.5,a,i[k]))-log(i[k])+ a*log(a/p[p>0])-(i[k]+a)*log(1+a/p[p>0])-log(1-(a/(a+p[p>0]))^a)))
  }
  prob
}

