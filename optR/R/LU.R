#' Solving system of equations using LU decomposition
#' 
#' @description The function solves Ax=b using LU decomposition (LUx=b). The function handles multple responses
#' @param A   : Input Matrix
#' @param b   : Response
#' @param tol : tolerance
#' @return U  : Upper part of the triangele is (U) and Lower part of the triangular is L (Diagnoal for the L matrix is 1) 
#' @return c  : Lc=b (where Ux=c)
#' @return beta : Estimates
#' @examples 
#' A<-matrix(c(0,-1,1, -1,2,-1,2,-1,0), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(0,0, 1), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="LU")
LU.optR<-function(A, b, tol=1e-7){
  method="LU"
  # Matrix check
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

    
  
  # LU decompose
  if(method=="choleski"){
    A<-choleskiDecomposition(A, tol)
    A<-A+t(A)
    diag(A)<-diag(A)/2
  } else
  {
    A<-LU.decompose(A, tol)
  }
  
  
  if(!is.null(b)){
    if(!is.matrix(b)) b<-as.matrix(b)
    nCOL<-ncol(b)
    y<-list()
    xEstimates<-list()
    for(i in 1:nCOL){
      # Estimate the y matrix
      y[[i]]<-forwardsubsitution.optR(A, b[, i])
      
      #Estimate coefficients of x solve: Ux=y
      xEstimates[[i]]<-optR.backsubsitution(A, y[[i]])
    }
    return(list("U"=A, "c"=y, "beta"=xEstimates, "seq"=seq))
  } else
  {
    return(list("U"=A, "seq"=seq))
  }
}

#' LU decomposition
#' 
#' @description The function decomposes matrix A into LU with L lower matrix and U as upper matrix
#' @param A   : Input Matrix
#' @param tol : tolerance
#' @return A  : Transformed matrix with Upper part of the triangele is (U) and Lower part of the triangular is L (Diagnoal for the L matrix is 1) 
#' @examples 
#' A<-matrix(c(0,-1,1, -1,2,-1,2,-1,0), nrow=3,ncol=3, byrow = TRUE)
#' Z<-optR(A, tol=1e-7, method="LU")
LU.decompose<-function(A, tol=1e-7) {
  nROW<-ncol(A)
  nCOL<-nROW
  
  # LU Decomposition
  for(i in 1:(nROW-1)){
    # Estimate the Lambda for
    lmbd<-matrix(unlist(lapply(1:nROW, FUN=optR.multiplyfactor, A, i)), nrow=nROW)
    
    # Minimize the iterations
    if(sum(abs(lmbd))!=0) {
      lmbdMatrix<-matrix(rep(lmbd[(i+1):nROW], each=nCOL), nrow = (nROW-i), byrow=T)
      A[(i+1):nROW, (i+1):nROW] <- A[(i+1):nROW,(i+1):nROW]-(lmbdMatrix*matrix(rep(A[i, 1:nROW], each=nROW-i), nrow = (nROW-i), byrow=F))[,(i+1):nROW]
      A[(i+1):nROW,i]=lmbd[(i+1):nROW]
    }
  }
  return(A)
}