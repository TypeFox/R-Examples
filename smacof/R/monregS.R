# secondary approach to ties

`monregS` <-
function(x, y, w = rep(1,length(x)), block = weighted.mean)
{
#x ... observed distance matrix
#y ... Guttman transformed distances
#w ... weights

  wag <- tapply(w,x,sum)
  yag <- tapply(y,x,mean)
  xag <- tapply(x,x,mean)
  o <- order(xag)
  r <- order(o)
  e <- pavasmacof(yag[o],wag[o])[r]                  #call pava
  return(ifelse(outer(x,xag,"=="),1,0)%*%e)
}

