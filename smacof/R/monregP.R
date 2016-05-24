# primary approach to ties
`monregP` <-
function(x, y, w = rep(1,length(x)), block = weighted.mean)
{
#x ... observed distances (diss)
#y ... Guttman transformed distances
#w ... weights

  o <- order(x,y)
  r <- order(o)
  return(pavasmacof(y[o],w[o])[r])
}

