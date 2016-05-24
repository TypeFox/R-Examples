mop.q = function (x, k, p, q, method=c("MOP", "RBMOP")) 
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
  
    
  if ( q < 0 || q > 1 || !is.numeric(q)) 
  {
    stop("q must lie between 0 and 1.")
  }
  
  method = match.arg(method)
  
  
  
  EVI.est = mop(x,k,p,method)
  EVI.est = EVI.est$EVI
  
  
  osx = sort(x) 
   dk = length(k)
   dp = length(p)
  
 mopVaR = matrix(NA,dk,dp)
  
    for(j in 1:dp)
    { 
      mopVaR[,j] = (osx[n-k]) *  (k/(n*q))^EVI.est[,j]
    }
  
  
  colnames(mopVaR) = p
  rownames(mopVaR) = k
  
  
  return(list(EVI=EVI.est,HQ=mopVaR)) 
   
}




  
