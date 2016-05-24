mop = function (x, k, p, method=c("MOP", "RBMOP")) 
{ 
  
  n = length(x) 
  
  # Checking for plausible inputes
  
  if (n  < 2) 
  {
    stop("Data vector x with at least two sample points required.")
  }
  
  if (is.null(p) || any(is.na(p)))
  {
    stop("p is not specified")
  }  
  
  
  if (is.null(k) || any(is.na(k)))
  {
    stop("k is not specified")
  }
  
  if (any(k < 1) || any(k >n) || any(k == n) || !is.numeric(k) || k != as.integer(k))
  {
    stop("Each k must be integer and greater than or equal to 1 and less than sample size.")
  }
  
  if (!is.numeric(x)) 
  {
    stop("Some of the data points are not real.")
  }
  
  if (!is.numeric(p)) 
  {
    stop("p must be a real.")
  }
  
  
  # Order statistics
  
  osx <- sort(x[x > 0])
  
  method = match.arg(method)
  
  
  if (method == "MOP") 
  {
    EVI.est = mop.MOP(osx,k,p)
  }
  
  else if (method == "RBMOP")
  {
    EVI.est = mop.RBMOP(osx,k,p)
  } 
  
  colnames(EVI.est$EVI) = p
  rownames(EVI.est$EVI) = k
  
  return(EVI.est)
  
}



# mean of order p extreme value index estimate

mop.MOP = function(osx,k,p)
{
  
  
  n = length(osx) 
  
  dk = length(k)
  dp = length(p)
  km = max(k)
  nk = km + 1
  
  est = matrix(NA,dk,dp)
  
  
    for(j in 1:dp)  
      
    # Hill estimation
    if (p[j] == 0)
    {
      
      losx = rev(log(osx))       
      losx = losx[1:nk]
      tmp  = cumsum(losx)/(1:nk)
      estc = tmp[-nk]-losx [-1]  
      est[,j]= estc[k]
      
    }
  
  else 
  {    
    tosx = rev(osx)
    tosx = (tosx[1:nk])^p[j]
    tmp = cumsum(tosx)/(1:nk)
    estc = tmp[-nk]/tosx[-1]
    estc = (1- estc^-1)/ p[j] 
    est[,j] = estc[k]
    
  } 
  return(list(EVI=est))
  
}  


# Reduced bias mean of order p extreme value index estimate

mop.RBMOP = function(osx,k,p)
{ 
  
  
  n = length(osx)
  losx = log(osx)       # log of order statistics
  dk = length(k)
  dp = length(p)
  
  est = matrix(NA,dk,dp)
  
  H =  mop.MOP(osx,k,p)
  H =  H$EVI  # EVI estimates without reduced bias
  
  # Second order parameteres  estimates
  rhoest  = mop.rho(osx)  
  betaest = mop.beta(losx,rhoest) 
  
  for(j in 1:dp)
  est[,j] = H[,j]* (1 -   (  (betaest* (1-p[j]*H[,j])) / (1-rhoest-p[j]*H[,j])  )* (n/k)^(rhoest)  )
    
  return(list(EVI=est,rho=rhoest,beta=betaest))
}  



# Rho estimation
mop.rho = function(osx)
{
  
  losx = log(osx) 
  n = length(losx)
  
  krho = floor(n^(0.995)):floor(n^(0.999))
  krhom = max(krho)
  nkrho = krhom + 1
  
  M = matrix(NA,length(krho),3)
  
  tau0 = numeric(length(krho))
  tau1 = numeric(length(krho))   
  W0 = numeric(length(krho))
  W1 = numeric(length(krho))
  
  losx = rev(losx)       
  losx = losx[1:nkrho]

  
  estc1 = mop.MOP(osx,krho,p=0)
  M[,1] = estc1$EVI
  
  c12 = cumsum(losx^2)/(1:nkrho)
  c22 = (losx)^2
  c32 = cumsum(losx)/(1:nkrho)
  estc2 = c12[-nkrho] + c22[-1] -2*c32[-nkrho]* losx[-1]
  M[,2]=estc2[krho]
  
    
  c13 = cumsum(losx^3)/(1:nkrho)
  c23 = (losx)^3
  estc3 = c13[-nkrho] - c23[-1] - 3* (losx[-1]*c12[-nkrho] - c32[-nkrho]*c22[-1])
  M[,3]=estc3[krho]
   
  
  
  W0 =  ( log(M[,1]) - (1/2)* log(M[,2]/2) )/ ( (1/2)* log(M[,2]/2) - (1/3)* log(M[,3]/6) )    
  W1 =  ( M[,1] - (M[,2]/2)^(1/2)  ) / ( (M[,2]/2)^(1/2) -  (M[,3]/6)^(1/3) )

 
  tau0 = -abs( 3*(W0-1)/(W0 -3) )
  tau1= -abs( 3*(W1-1)/(W1 -3) )    
 
 
  tau0median = median(tau0)
  tau1median = median(tau1)
  
  
  sumtau0 = as.numeric((tau0 - tau0median) %*% (tau0 - tau0median))
  sumtau1 = as.numeric((tau1 - tau1median) %*% (tau1 - tau1median))
  
    
  if((sumtau0 < sumtau1) || (sumtau0 == sumtau1) )
  {
    return(tau0[length(krho)])   
    
  }
  
  else
  {
    return(tau1[length(krho)])
  }
  
  
}

# Beta estimation
mop.beta = function(losx,rhoest)
{
  n = length(losx)  
  k1 = floor(n^(0.999))
  
  v = c(0,rhoest,2*rhoest)
  D = numeric(length(v))
  
  c1 = ((1:k1)/k1)^(-v[1])
  c2 = ((1:k1)/k1)^(-v[2])
  c3 = ((1:k1)/k1)^(-v[3])
  
 
  
  c =(1:k1)* ( (losx[n:(n-k1+1)]) - (losx[(n-1):(n-k1)]) )
   
  D[1] = sum(c1*c)/k1
  D[2] = sum(c2*c)/k1
  D[3] = sum(c3*c)/k1
  
  
  d = sum(c2)/k1
  
  
  betaest  =  (k1/n)^(v[2])*(d*D[1] -D[2])/(d*D[2]-D[3]) 

  return(betaest)
  
  
}
