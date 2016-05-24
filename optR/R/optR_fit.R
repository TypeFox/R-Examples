#' Fitter function for Linear/Non-linear system with form Ax=b  
#' 
#' optR.fit is fit function for determing x for System with form Ax=b
#' @param x     : Input matrix
#' @param y     : Response is matrix
#' @param method  : "gauss" for gaussian elimination and "LU" for LU factorization
#' @param iter  : Number of Iterations
#' @param tol   : Convergence tolerance
#' @param ...   : S3 Class
#' @return U    : Decomposed matrix for Gauss-ELimination Ax=b is converted into Ux=c where U is upper triangular matrix for LU decomposition U contain the values for L & U decomposition LUx=b   
#' @return c    : transformed b & for LU transformation c is y from equation Ux=y
#' @return estimates  : Return x values for linear system
#' @return seq        : sequence of A matrix re-ordered
#' @examples
#' # Solving equation Ax=b
#' A<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(-14,36, 6), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gauss") 
#' 
#' # Solve Linear model using LU decomposition (Supports Multi-response)
#' Z<-optR(A, b, method="LU")
#' 
#' # Solving the function using numerical method
#' Z<-optR(A, b, method="cgm")
optR.fit<-function(x, y=NULL, method=c("gauss, LU, gaussseidel", "cgm"), iter=500, tol=1e-7, ...){
  optR<-list()
  if(method=="gauss" & nrow(x)==ncol(x)){
    Z<-optR.gauss(x, y, tol)
    optR$U<-Z$U
    optR$c<-Z$c
    optR$beta<-Z$beta
  } else if (method=="LU" & nrow(x)==ncol(x)) {
    Z<-LU.optR(x, y, tol)
    optR<-assign.values(Z)
  } else if (method=="gaussseidel" & nrow(x)==ncol(x)) {
    Z<-gaussSeidel(x, y, iter=iter, tol=tol)
    optR$beta<-Z$optimal
    optR$initial.values<-Z$initial
    optR$relaxationFactor<-Z$relaxationFactor
    optR$itr.conv<-Z$itr.conv
  } else if (method=="cgm" & nrow(x)==ncol(x)) {
    Z<-cgm(x, y, iter=iter, tol=tol)
    optR$beta<-Z$optimal
    optR$initial.values<-Z$initial
    optR$itr.conv<-Z$itr.conv
    optR$conv<-Z$conv
  } else if(method=="choleski" & nrow(x)==ncol(x)) {
    Z<-choleskilm(x, y, tol=tol)
    optR<-assign.values(Z)
  } else
  {
    stop("Method does found...")
  }
  class(optR)<-"optR"
  optR
}


# Function assigns values for LU and Choleski decomposition
assign.values<-function(Z){
  optR<-list()
  optR$U<-Z$U
  optR$c<-Z$c
  optR$beta<-Z$beta
  optR$seq<-Z$seq
  return(optR)
}
