# ' Compute the log-likelihood of a gaussian mixture value due to missing values in a regression
# ' @param X is the vector of the regressors
# ' @param Y is the value of the response variable (scalar)
# ' @param B is the vector of regression parameters, including the intercept in first position  (zero if no intercept)
# ' @param sigma is the standard deviation of the residual of the subregression
# ' @param M is the indicatrice vector (submatrix) of the missing values
# ' @param mixmod is result of calcul_BIC_mixmod2.0(X=dataset,nbclustmax=nbclustmax,details=T)
# ' @param log boolean to define if you want the likelihood or the log-likelihood
# ' @param intercept boolean to define if there is an intercept in the model (default=T)
GM_Loglikelihood<-function(Y=Y,X=X,B=B,sigma=sigma,M=NULL,mixmod,log=T,intercept=T){
  if(is.null(M)){
    M=0*X
    M[is.na(M)]=1
    X[is.na(X)]=0
  }
  return(.Call( "GM_likelihood",Y,X,B,sigma,M,mixmod$nbclust,mixmod$details,log,intercept, PACKAGE = "CorReg"))
}


