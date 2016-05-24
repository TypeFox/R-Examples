mop.AREFF = function (x, k, p)   
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
  
  if (!is.numeric(x)) {
    stop("Some of the data points are not real.")
  }
  
  if (!is.numeric(p)) 
  {
    stop("p must be a real.")
  }
  
   
  
  EVI.est = mop(x,k,p,"MOP")
  EVI.est = EVI.est$EVI
  
  RBEVI   = mop(x,k,p,"RBMOP")
    
  rhoest = RBEVI$rho

 if(any(p*EVI.est>0.5))
 {
    stop("Condition on EVI and p  is not satified.")
 }
 
  dk = length(k)
  dp = length(p)
  
  c1 = numeric(dk)
  c2 = numeric(dk)
  
  AREFF = matrix(NA,dk,dp)
  
  for(j in 1:dp)
  {  
         c1 = ( sqrt(1-2*p[j]*EVI.est[,j])/(1-p[j]*EVI.est[,j])) ^(-2*rhoest)   
         c2 = abs( (1-p[j]*EVI.est[,j]-rhoest)/((1-rhoest)*(1-p[j]*rhoest)) )
   AREFF[,j]= (c1*c2)^(1/(1-2*rhoest))
  }  
  


  colnames(AREFF) = p
  rownames(AREFF) = k
    
  return(list(AREFF=AREFF))  
  
}