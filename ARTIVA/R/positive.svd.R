positive.svd <-
function(m, tol)
{
  s = svd(m)
  
  if( missing(tol) ) 
      tol = max(dim(m))*max(s$d)*.Machine$double.eps
  Positive = s$d > tol

  return(list(
      d=s$d[Positive],
      u=s$u[, Positive, drop=FALSE],
      v=s$v[, Positive, drop=FALSE]
      ))
}
