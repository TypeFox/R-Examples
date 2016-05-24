fast.svd <-
function(m, tol)
{  
  n = dim(m)[1]
  p = dim(m)[2]
 
 
  EDGE.RATIO = 2 # use standard SVD if matrix almost square
  if (n > EDGE.RATIO*p)
  {
     return(psmall.svd(m,tol))
  }
  else if (EDGE.RATIO*n < p)
  {  
     return(nsmall.svd(m,tol)) 
  }
  else # if p and n are approximately the same
  {
     return(positive.svd(m, tol))
  }
}
