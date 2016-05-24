`centerX` <- function(x,w)
{
  apply(x,2,function(z) z-weighted.mean(z,w))
}
