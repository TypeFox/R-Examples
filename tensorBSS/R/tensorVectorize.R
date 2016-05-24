tensorVectorize <-
function(x){
  rplusone <- length(dim(x))
  apply(x, rplusone, c)
}
