other.q = function (x,k,q, method=c("MO", "GH","MM"))   
{ 
  
  n = length(x) 
  
  # Checking for plausible inputes
  
  if (n  < 2) {
    stop("Data vector x with at least two sample points required.")
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
    stop("some of the data points are not real.")
  }
  
    
  if ( q < 0 || q > 1 || !is.numeric(q)) 
  {
    stop("q must lie between 0 and 1.")
  }
  
  method = match.arg(method)
  
  
  EVI.est = other.EVI(x,k,method)
  EVI.est = EVI.est$EVI
  
  osx = sort(x) 
  
  dk = length(k)
  otherVaR = numeric(dk)  
  
  otherVaR = (osx[n-k]) *  (k/(n*q))^EVI.est
  
  otherVaR = as.matrix(otherVaR)
  rownames(otherVaR)=k 
  return(list(EVI=EVI.est, HQ=otherVaR))  
    
}
