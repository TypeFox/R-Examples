# ' calcul du BIC theta
# ' 
BicTheta<-function(X=X,Y=Y,beta=beta,intercept=T){
  res=.Call( "BicTheta",matrix(as.double(X),ncol=ncol(X)),as.double(Y),beta,intercept, PACKAGE = "CorReg")
  return(res$BIC)
}