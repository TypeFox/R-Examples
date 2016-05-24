# ' Cross validation
# ' @export
# ' @param X covariates matrix (double)
# ' @param Y response variable
# ' @param K number of classes
# ' @param intercept (boolean) with or without an intercept
# ' @param methode the methode used by OLS.
# ' @param groupe a vector to define the groups used for cross-validation (to obtain a reproductible result)

CVMSE<-function(X=X,Y=Y,K=K,intercept=T,methode=1,groupe=NULL){
  K=abs(K)
  K=min(K,nrow(X))
  if(is.null(groupe)){
    groupe=rep(0:(K-1),length.out=nrow(as.matrix(X)))
    groupe=sample(groupe)
  }
  res=.Call( "CVMSE",X,as.double(Y),K,intercept,methode,groupe, PACKAGE = "CorReg")
  return(res$MSE)
}