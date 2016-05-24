#' Linear combinations of submatrices of an array
#' 
#' Computes a matrix of expected values based on an array X of predictors and a
#' vector beta of regression coefficients.
#' 
#' 
#' @usage Xbeta(X, beta)
#' @param X an n by n by p array
#' @param beta a p by 1 vector
#' @return An n by n matrix
#' @author Peter Hoff
#' @export Xbeta
Xbeta <-
function(X,beta)
{
  XB<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2] )
  for(k in seq(1,length(beta),length=length(beta))){XB<-XB + beta[k]*X[,,k]}
  XB
}
