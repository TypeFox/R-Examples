Px <-
function(X)
{
  if(is.vector(X))
  {
    X=as.matrix(X)
  }
  
  res=X%*%solve(t(X)%*%X)%*%t(X)
  return(res)
}
