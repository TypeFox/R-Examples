#' simple MSE function
#' @description
#' This function computes the MSE (Mean Squared Error) of prediction associated to a vector of coefficients \code{A} used to predict a response variable \code{Y} by linear regression on \code{X}, with an intercept or not.
#' @export
#' @param Y the response variable (vector)
#' @param X the dataset (matrix of covariates)
#' @param A the vector of coefficients
#' @param intercept (boolean) to add a column of 1 to \code{X} if \code{A} contains an intercept and \code{X} doesn't.
#' @return the Mean Squared Error observed on \code{X} when using \code{A} coefficients to predict \code{Y}.
#' 
#'  @examples
#'  require(CorReg)
#'    #dataset generation
#'    base=mixture_generator(n=15,p=5,valid=100,scale=TRUE)
#'    X_appr=base$X_appr #learning sample
#'    Y_appr=base$Y_appr#response variable
#'    X_test=base$X_test#validation sample
#'    Y_test=base$Y_test#response variable (validation sample)
#'    A=lm(Y_appr~X_appr)$coefficients
#' MSE_loc(Y=Y_appr,X=X_appr,A=A)#MSE on the learning dataset
#' MSE_loc(Y=Y_test,X=X_test,A=A)#MSE on the validation sample
MSE_loc<-function(Y=Y,X=X,A=A,intercept=T){
  if(intercept){
    X=as.matrix(cbind(1,X))
  }
  return(mean((Y-X%*%A)^2))
}