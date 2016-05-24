# ' clean Z (if BIC improved)
# ' @export
# '@param X the dataset
# '@param Z binary adjacency matrix of the structure (size p)
# ' @param Bic_null_vect the BIC of the null hypothesis (used for independent variables)
# ' @param methode parameter for OLS (matrix inversion) methode_BIC  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
# ' @param star boolean defining wether classical BIC or BIC* is computed
# '@param plot if TRUE returns the vector of BIC for each step
# '@param verbose 0:none, 1:BIC,step and complexity when best BIC found 2:BIC, step, complexity, nb candidates and best candidate when best BIC found
# '
cleanZ<-function(X=X,Z=Z,Bic_null_vect=Bic_null_vect,methode=1,plot=F,verbose=1,star=FALSE){
   X=1*as.matrix(X)
   res=.Call( "cleancolZ",X,Z,Bic_null_vect,methode,plot,verbose,star, PACKAGE = "CorReg")
  res=.Call( "cleanZ",X,res$Z,Bic_null_vect,methode,plot,verbose,star, PACKAGE = "CorReg")
  return(res)
}