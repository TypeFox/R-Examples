#' Function fits linear model using Choleski Decomposition
#' 
#' @description Function fits a linear model using Choleski Decomposition for positive definate matrix
#' @param A : Input matrix
#' @param b : Response matrix
#' @param   tol : Tolerance
#' @return U  : Upper part of the triangele is (U) and Lower part of the triangular is L (Diagnoal for the L matrix is 1) 
#' @return c  : Lc=b (where Ux=c)
#' @return beta : Estimates
#' 
#' examples
#' A<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(-14,36, 6), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="choleski") # Solve Linear model using Gauss Elimination
choleskilm<-function(A, b, tol=1e-7){
  if(!is.matrix(b) & !is.null(b)) b<-as.matrix(b)
  
  # Matrix re-ordering to improve the factors
  A_ordered<-opt.matrix.reorder(A, tol)
  A<-A_ordered$A
  seq<-A_ordered$b.order
  rm(A_ordered)
  
  # Re-order b if it is not NULL
  if(!is.null(b)){
    b<-b[seq, ]  
  }
  # choleski decomposition
  A<-choleskiDecomposition(A, tol)
  
  if(!is.null(b)){
    if(!is.matrix(b)) b<-as.matrix(b)
    nCOL<-ncol(b)
    y<-list()
    xEstimates<-list()
    for(i in 1:nCOL) {
      
      # Estimate the y matrix
      y[[i]]<-forwardsubsitution.Cholesky(A, b[, i])
      
      #Estimate coefficients of x solve: Ux=y
      xEstimates[[i]]<-optR.backsubsitution(t(A), y[[i]])
    }    
    return(list("U"=A, "c"=y, "beta"=xEstimates, "seq"=seq))
    } else
  {
    return(list("U"=A, "seq"=seq))
  }
}


#' Function for Choleski Decomposition
#' 
#' @description Function perform choleski decomposition for positive definate matrix (A=LL')
#' @param A :Input Matrix
#' @param tol : Tolerance
#' @return L: Decomposition matrix
#' @export
#' @examples
#' A<-matrix(c(4,-2,2, -2,2,-4,2,-4,11), nrow=3,ncol=3, byrow = TRUE)
choleskiDecomposition<-function(A, tol=1e-7){
  nROW<-ncol(A)
  L<-matrix(rep(0, each=nROW*nROW), nrow = nROW, byrow=T)
  for(i in 1:nROW){
    
    # Fill Lower triangular matrix for diagnoal matrix
    Aii<-A[i,i]-sum(L[i, 1:i]*L[i, 1:i])
    if(Aii<0){
      stop("Matrix no positive definate")
    } else 
    {
      L[i,i]<-sqrt(Aii)  
    }
    
    # Fill Lower triangular matrix for non-diagnoal matrix
    if((i+1)<=nROW) {
      for(k in (i+1):nROW){
        L[k,i]<-(A[k,i]-sum(L[k,1:i]*L[i,1:i]))/L[i,i]
      }  
    }
  }
  return(L)
}
