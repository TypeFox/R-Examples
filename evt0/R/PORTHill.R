PORT.Hill = function (x, k, q, method=c("PMOP", "PRBMOP")) 
  
{ 
    
  n = length(x) 
  
  nq = floor(n*q) + 1
  
  # Checking for plausible inputes
  
  if (nq < 2 ) 
  {
    stop("Insufficient data vector, check sample size and q.")
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
  
  
  if ( q < 0 || q > 1 || !is.numeric(q)) 
  {
    stop("q must lie between 0 and 1.")
  }
  
  
  method = match.arg(method)
  
  
  if (method == "PMOP") 
  {
    PORTEVI = PORT.MOP(x,k,q)
  }
  
  else if (method == "PRBMOP")
  {
    PORTEVI = PORT.RBMOP(x,k,q)
  }  
  
 # PORTEVI= as.matrix(PORTEVI$PORT.EVI)
  #rownames(PORTEVI)=k 
  return(PORTEVI)  
  
   
}

# PORT Hill estimator

PORT.MOP = function(x,k,q)
{
  
  n = length(x)
  osx = sort(x) 
  nq = floor(n*q) + 1
  
  dk = length(k)
  est = numeric(dk) 
  
  km = max(k)
  nk = km + 1
  
  rosx  = rev(osx)  
  rosxq = rosx[1:nk] -osx[nq]
  
  losx  = log(rosxq)
  
  tmp   = cumsum(losx)/(1:nk)
  estc  = tmp[-nk]-losx [-1]  
  est   = estc[k]
      
     
    
  return(list(PORT.EVI=est))
    
}


# quasi-PORT Hill estimator

PORT.RBMOP = function(x,k,q)
{ 
  n = length(x)
  dk = length(k)
 est = numeric(dk) 
 
  
  H =  PORT.MOP(x,k,q)
  H =  H$PORT.EVI
   
  S =  mop(x,k[1],p=0,"RBMOP")
  S1= S$rho
  S2= S$beta
   
  
  est = H* (1 -   ( S2 / (1-S1)  )* (n/k)^(S1)    )
  
  return(list(PORT.EVI =est, rho =S1, beta=S2))
}  
