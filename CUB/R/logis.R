#' @title The logistic transform
#' @description Create a matrix YY binding array Y with a vector of ones, placed as the first column of YY. 
#' It applies the logistic transform componentwise to the standard matrix multiplication between YY and param.
#' @aliases logis
#' @usage logis(Y, param)
#' @export logis
#' @param Y A generic matrix or one dimensional array
#' @param param Vector of coefficients, whose length is NCOL(Y) + 1 (to consider also an intercept term)
#' @return Return a vector whose length is NROW(Y) and whose i-th component is the logistic function
#' at the scalar product between the i-th row of YY and the vector param.
#' @keywords utilities
#' @examples
#' n<-50 
#' Y<-sample(c(1,2,3),n,replace=TRUE) 
#' param<-c(0.2,0.7)
#' logis(Y,param)

logis <-
function(Y,param){
  YY<-cbind(1,Y)          # add 1's first column to matrix Y
  return(1/(1+exp(-YY%*%param)))
}
