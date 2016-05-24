wmean2 <-
function(zmu,y) {
  l <- length(zmu)
  weighted.mean( (y-zmu[l])^2, zmu[1:(l-1)] )
}
