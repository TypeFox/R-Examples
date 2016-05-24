# ' Bic of one regression
# ' 

BicLoc<-function(X=X,Y=Y,intercept=T,sigma=T,SumSquare=F,methode=1){
  ret= .Call( "BicLoc", X,Y,intercept,sigma,SumSquare,methode, PACKAGE = "CorReg")
  return(ret)
}