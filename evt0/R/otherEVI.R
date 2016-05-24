other.EVI = function (x, k, method=c("MO", "GH","MM")) 
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
  
 
  
  
  method = match.arg(method)
  
  
  if (method == "MO") 
  {
    EVI.est = mo(x,k)
  }
  
  else if (method == "GH")
  {
    EVI.est = gh(x,k)
  }  
  
  else if (method == "MM")
  {
    EVI.est = mm(x,k)
  }  
  
  EVI.est = as.matrix(EVI.est)
  rownames(EVI.est) = k 
   
  return(list(EVI=EVI.est))
  
  
}


# Moment(MO) estimator
mo = function(x,k)
{
   osx = sort(x) 
     n = length(osx)
  losx = log(osx)
    dk = length(k)
moment = numeric(dk)
     M = matrix(NA,dk,2)
  
     km = max(k)
     nk = km + 1
      
   losx = rev(losx)       
   losx = losx[1:nk]
   
     estc1 = mop(x,k,p=0,"MOP")
     M[,1] = estc1$EVI
     
               
     c12 = cumsum(losx^2)/(1:nk)
     c22 = (losx)^2
     c32 = cumsum(losx)/(1:nk)
   estc2 = c12[-nk] + c22[-1] - 2*c32[-nk]* losx[-1]
   M[,2] = estc2[k]
      
  moment = M[,1] + (1/2)*( 1 -  (  ( M[,2]/(M[,1])^2) - 1)^(-1) )
   #moment= M[,1] + 1 - (1/2)*( 1 -(M[,1])^2/ M[,2]) ^(-1)
     
   
   return(moment)
  
}  


# Generalised Hill (GH) estimator
gh = function(x,k)
{
  osx = sort(x) 
  n = length(osx)
  
  dk = length(k)
  ghest= numeric(dk)
  
  km = max(k) 
  
      c =  mop(x,1:km,p=0,"MOP")
      c = c$EVI
  
    evi = c[k]  
   levi = log(evi)  
   
   temp = cumsum(log(c))/(1:km)
  ghest = evi + (temp[k] - levi) 
   
  return(ghest)
  
} 

# Mixed Moment(MM) estimator
mm = function(x,k)
{
  osx = sort(x) 
  n = length(osx)
  km = max(k)
  nk = km + 1
  
  
  dk = length(k)
  mmest= numeric(dk)
  
  estc1 = mop(x,k,p=0)
  M     = estc1$EVI
    
  tosxv = rev(osx)
  tosxv = tosxv[1:nk]
   tosx = (tosxv)^(-1)
   tmp  = cumsum(tosx)/(1:nk)
   estc = tosxv[-1]*tmp[-nk]
      L = estc[k]
    
    phi = (M -L)/(L^2)    
  mmest = (phi-1)/(1 + 2*min(phi-1,0))
    
   return(mmest)
  
}   

