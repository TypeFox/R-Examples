#' Optimization & estimation based on Conjugate Gradient Method 
#' 
#' @description Function utilizes the Conjugate Gradient Method for optimization to solve equation Ax=b
#' @param A     : Input matrix
#' @param b     : Response vector
#' @param x     : Initial solutions
#' @param iter  : Number of Iterations
#' @param tol   : Convergence tolerance
#' @return optimal  : Optimal solutions
#' @return initial  : initial solution
#' @return itr.conv  : Number of iterations for convergence
#' @return conv  : Convergence array
#' @examples
#' A<-matrix(c(4,-1,1, -1,4,-2,1,-2,4), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(12,-1, 5), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="cgm", iter=500, tol=1e-7)
cgm<-function(A, b, x=NULL, iter=500, tol=1e-7){
  nROW<-nrow(A)
  if(is.null(x)) x<-matrix(rep(0, each=nROW), nrow = nROW, byrow=T)
  conv<-NULL
  xini<-x
  SearchDirection<-0
  Beta<-0
  for(i in 1:iter){
    xold<-x
    residual<-b-A%*%x
    if(i>1){
      Beta<--1L*(t(residual)%*%A%*%SearchDirection)/(t(SearchDirection)%*%A%*%SearchDirection)
    }
    SearchDirection<-residual + matrix(rep(Beta, each=nROW), nrow = nROW, byrow=T)*SearchDirection
    
    # CHECK if solution is optimized
    alphaParameter<-(t(SearchDirection)%*%residual)
    if(alphaParameter==0 & i>1){
      return(list("x"=x, "xini"=xini, "itr.conv"=i, "conv"=conv))
    }
    
    # Update alpha & x parameter
    alpha=alphaParameter/(t(SearchDirection)%*%A%*%SearchDirection)
    x<-x+matrix(rep(alpha, each=nROW), nrow = nROW, byrow=T)*SearchDirection
    dx<-sqrt(t(x-xold)%*%(x-xold))
    conv<-c(conv, dx)
    if(dx<=tol){
      return(list("optimal"=x, "initial"=xini, "itr.conv"=i, "conv"=conv))
    }
  }
  print("Optimization Failed to Converge...")
  return(list("x"=x, "xini"=xini, "itr.conv"=i, "conv"=conv))
}