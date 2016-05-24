# ' clean Z columns (if BIC improved)
# '@param X the dataset
# '@param Z binary adjacency matrix of the structure (size p)
# '@param Bic_null_vect vector of the BIC for each covariate
# '@param methode_BIC  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
# '@param star if TRUE BIC* is used instead of bic
# '@param plot if TRUE returns the vector of BIC for each step
# '@param verbose 0:none, 1:BIC,step and complexity when best BIC found 2:BIC, step, complexity, nb candidates and best candidate when best BIC found
# '@export
cleancolZ<-function(X=X,Z=Z,Bic_null_vect=Bic_null_vect,methode_BIC=1,plot=F,verbose=1,star=FALSE){
   X=1*as.matrix(X)
   res=.Call( "cleancolZ",X,Z,Bic_null_vect,methode_BIC,plot,verbose,star, PACKAGE = "CorReg")
  return(res)
}