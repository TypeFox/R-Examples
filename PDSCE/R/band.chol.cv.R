band.chol.cv <-
function(x, k.vec=NULL, method= c("fast", "safe"), nsplits=10, n.tr=NULL, quiet = TRUE)
{
  method=match.arg(method)
  n=dim(x)[1]
  p=dim(x)[2]
  if(is.null(n.tr))
    n.tr=round(n*(1-1/(log(n))))
  n.va=n-n.tr
  
  if(is.null(k.vec))
  {
    k.vec=0:min(c(n-2, p-1))
  }
  
  cv.loss = array(0, c(length(k.vec), nsplits))
  for ( ns in 1:nsplits) 
  {
    ind = sample(n)
    ind.tr = ind[1:n.tr]
    ind.va = ind[(n.tr+1):n]
    x.tr=scale( x[ind.tr,,drop=FALSE], scale=FALSE)
    s.va = cov(x[ind.va,,drop=FALSE]) * (n.va -1)/n.va
    
    for( i in 1:length(k.vec))
    {
      sigma=band.chol(x=x.tr,k=k.vec[i], centered=TRUE, method=method)
      cv.loss[i,ns] = sum ( (sigma - s.va)^2 )    
      if(!quiet) cat("Finished k =", k.vec[i], "in split", ns, "\n") 
    }
    if(!quiet) cat("Finished split", ns, "\n")       
  }  
  cv.err=apply(cv.loss, 1, sum)
  best.k = k.vec[which.min(cv.err)]
  sigma=band.chol(x=x,k=best.k, method=method)   

  return(list(sigma=sigma, best.k=best.k, cv.err=cv.err, k.vec=k.vec, 
              n.tr=n.tr)) 
}

