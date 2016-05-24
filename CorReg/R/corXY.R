# ' Correlations between X and Y
corXY<-function(X=X,Y=Y){
  p=ncol(X)
  return(as.numeric(cor(cbind(X,Y))[-(p+1),p+1]))
}