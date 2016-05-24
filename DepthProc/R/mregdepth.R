mregdepth <- function(coef, x, y, ndir = 1000)
{
  if(is.null(dim(x))) x = matrix(x, ncol = 1)

  yest = cbind(1,x)%*%coef
  res = y - yest
  u = runifsphere(ndir,length(coef))
  rdtmp = as.vector(res)/t(u%*%t(cbind(1,x)))<0
  min(colSums(rdtmp))
}
