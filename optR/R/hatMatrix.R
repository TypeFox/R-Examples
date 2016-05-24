#' Function determines the Hat matrix or projection matrix for given X
#' 
#' @description Function hatMatrix determines the projection matrix for X from the form yhat=Hy. The projection matrix defines the influce of each variable on fitted value
#' The diagonal elements of the projection matrix are the leverages or influence each sample has on the fitted value for that same observation. 
#' The projection matrix is evaluated with I.I.D assumtion ~N(0, 1)
#' @param X : Input Matrix
#' @param covmat : covariance matrix for error, if the error are correlated for I.I.D covmat will be NULL matrix
#' @return X: Projection Matrix
#' @export
#' @examples 
#' X<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' covmat <- matrix(rnorm(3 * 3), 3, 3)
#' H<-hatMatrix(X)
#' H<-hatMatrix(X, covmat)
#' diag(H)
hatMatrix<-function(X, covmat=NULL){
  if(!is.matrix(X)) X<-as.matrix(X)
  if(!is.null(covmat)){
    if(!is.matrix(covmat)) covmat<-as.matrix(covmat)
    if(dim(covmat)[1]!=dim(covmat)[2]) {
      stop("Covariance matrix no square")
    }
    A<-t(X)%*%inv.optR(covmat)%*%X
    A<-inv.optR(A)
    X<-X%*%A%*%t(X)%*%inv.optR(covmat)
  } else
  {
    A<-t(X)%*%X
    A<-inv.optR(A)
    X<-X%*%A%*%t(X)
  }
  return(X)
}