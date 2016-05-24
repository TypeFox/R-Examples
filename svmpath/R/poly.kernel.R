"poly.kernel" <-
  function(x, y=x, param.kernel = 1,...)
{
  if(is.null(param.kernel))
    param.kernel <- 1
  K= if(param.kernel == 1)
    x %*% t(y)
  else (x %*% t(y) + 1)^param.kernel
 if(param.kernel==1)  attr(K,"linear")=TRUE
  K
}
