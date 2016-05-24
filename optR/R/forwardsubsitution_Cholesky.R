forwardsubsitution.Cholesky<-function(L, b) {  
  nROW<-nrow(L)
  # Forward subsitution
  y<-b
  if(L[1,1]!=0) y[1]<-b[1]/L[1,1]
  for(k in 2:nROW) {
    if(L[1,1]!=0) {
      y[k]<-(y[k]-sum(L[k,1:(k-1)]*y[1:(k-1)]))/L[k,k]  
    } else
    {
      y[k]<-(y[k]-sum(L[k,1:(k-1)]*y[1:(k-1)])) # NAN values replacement
    }
    
  }
  return(y)
}