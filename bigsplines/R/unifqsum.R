unifqsum <-
function(x){
  x <- matrix(x,ncol=2)
  c(x[1,1],sum(x[,2]),length(x[,1]))
}
