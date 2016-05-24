otherPORT.q = function (x, k, q1, q2, method=c("MO", "GH","MM"))   
{ 
    
  n = length(x) 
  
  nq = floor(n*q1) + 1
  
  # Checking for plausible inputes
  
  if (nq < 2 ) 
  {
    stop("Insufficient data vector, check sample size and q1.")
  }
  
  if (is.null(k) || any(is.na(k)))
  {
    stop("k is not specified")
  }
  
  if (any(k < 1) || any(k > n-nq -1) || any(k == n-nq -1)  || !is.numeric(k)   || k != as.integer(k)  )
  {
    stop("k must be greater than or equal to 1 and less than exceedance sample size.")
  }
  
  if (!is.numeric(x)) {
    stop("some of the data points are not real.")
  }
  
  
  if ( q1 < 0 || q1 > 1 || !is.numeric(q1)) 
  {
    stop("q1 must lie between 0 and 1.")
  }
  
  if ( q2 < 0 || q2 > 1 || !is.numeric(q2)) 
  {
    stop("q2 must lie between 0 and 1.")
  }
  
  method = match.arg(method)
  
  osx = sort(x) 
  
  EVI.est = other.EVI(x,k,method)
  EVI.est = EVI.est$EVI
  
  dk = length(k)
  VaRest = numeric(dk)
  
 
  VaRest =   (osx[n-k]- osx[nq] )  * (k/(n*q2))^EVI.est + osx[nq]
  
  VaRest = as.matrix(VaRest)
  rownames(VaRest) = k 
  return(list(EVI=EVI.est, HQ=VaRest))  
  
}