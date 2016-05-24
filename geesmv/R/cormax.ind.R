cormax.ind <-
function(n){
  matrix<-matrix(0,nrow=n,ncol=n)
  diag(matrix)<-rep(1,n)
  return(matrix)
}
