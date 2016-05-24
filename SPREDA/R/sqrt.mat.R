sqrt.mat <-
function(x)
{
  aa=eigen(x)
  ee=aa[[1]]
  res=aa[[2]]%*%diag(sqrt(ee))%*%t(aa[[2]])
  return(res)
}
