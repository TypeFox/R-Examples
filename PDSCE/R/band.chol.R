band.chol <-
function(x,k, centered=FALSE, method=c("fast", "safe"))
{
  n=dim(x)[1]
  method=match.arg(method)  
  if (!centered)
  {
    mx=apply(x, 2, mean)
    x=scale(x, center=mx, scale=FALSE)
  }
  if(method=="fast")
  {
    sigma=fast.band(x=x, k=k)
  } else
  {
    sigma=safe.band(x=x, k=k)
  }
  return(sigma)
}

