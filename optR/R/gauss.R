#' gauss to solve linear systems 
#' 
#' Function solves linear systems using Gauss Elimination. The function solves equation of form Ax=b to Ux=c (where U is upper triangular matrix)
#' @param A       : Input Matrix
#' @param b       : Response
#' @param method  : To be used to perform factorization 
#' @param tol     : Tolerance
#' @return U      : Upper triangular matrix
#' @return c      : Transformed b
#' @return beta   : Estimates
#' @examples 
#' A<-matrix(c(0,-1,1, -1,2,-1,2,-1,0), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(0,0, 1), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gauss")
optR.gauss<-function(A, b, tol=1e-7){
  if(!is.matrix(A)) A<-as.matrix(A)
  if(!is.matrix(b)) b<-as.matrix(b)
  
  
  nCOL<-ncol(A)
  nROW<-nrow(A)
  
  # Matrix re-ordering to improve the factors
  if(prod(diag(A))==0){
    A_ordered<-opt.matrix.reorder(A, tol)
    A<-A_ordered$A
    b<-b[A_ordered$b.order, ]
    rm(A_ordered)
  }
  
  for(i in 1:(nROW-1)){
    
    # Estimate the Lambda
    lmbd<-matrix(unlist(lapply(1:nROW, FUN=optR.multiplyfactor, A, i)), nrow=nROW)
    
    # Raise Error if Inf Lambda identified
    if(sum(abs(lmbd))==Inf){
      stop("Pivot values wrongly ordered: Retuned Inf Value for Lambda")
    }
    
    # Minimize the iterations
    if(sum(abs(lmbd))!=0 ) {
      lmbdMatrix<-matrix(rep(lmbd[(i+1):nROW], each=nCOL), nrow = (nROW-i), byrow=T)
      A[(i+1):nROW,] <- A[(i+1):nROW,]-lmbdMatrix*matrix(rep(A[i, 1:nROW], each=nROW-i), nrow = (nROW-i), byrow=F)
      b[(i+1):nROW] <- b[(i+1):nROW] - lmbd[(i+1):nROW]*b[i]
    } 
  }
  
  # value less than machine precision assign zero
  A<-machinePrecision(A) # Machine precision error
  
  #Estimate coefficients
  xEstimates<-optR.backsubsitution(A, b)
  
  return(list("U"=A, "c"=b, "beta"=xEstimates))
}



#' Function to Re-order the matrix to make dominant diagnals
#' 
#' @description Function re-order the matrix to make matrix pivot.diag at each iteration
#' @param A         : Input Matrix
#' @param tol       : tolerance
#' @return A        : Updated Matrix
#' @return b.order  : Order sequence of A updated  matrix
#' @examples
#' A<-matrix(c(0,-1,1, -1,2,-1,2,-1,0), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(0,0, 1), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gauss")
opt.matrix.reorder<-function(A, tol=1e-16){
  # Initialize
  nROW<-nrow(A)
  b.order<-seq(1, nROW, by=1)
  
  # Assign max absolute value for column (Scale fator: absolute max  value of each rows)
  S<-apply(A, 1, FUN=function(x) max(abs(x)))
  for(k in 1:(nROW-1)){
    
    # find the row containint the largest relative size
    p<-which.max(abs(A[k:nROW, k]) /S[k:nROW])+k-1
    if(A[p, k]<tol){ 
      warning('Singular Matrix')
    }
    
    # Swapping rows
    if(p!=k){
      tmp<-A[k, ]
      A[k, ]<-A[p, ]
      A[p, ]<-tmp
      
      # Sequence of change
      tmp<-b.order[k]
      b.order[k]<-b.order[p]
      b.order[p]<-tmp
    }
  }
  #print(b.order)
  return(list("A"=A, "b.order"=b.order))
}



#' Function to estimate lambda
#' 
#' Function esimates the lambda or multiplier factor for Elimination using the pivot row/column
#' @param rowindex    : Row Index for the row to be used
#' @param pivotindex  : Column index for the pivot
#' @param A           : Input matrix
#' @return lambda     : Lambda
optR.multiplyfactor<-function(rowindex, A, pivotindex) {
  pivotVal<-A[rowindex,pivotindex]
  if(pivotVal!=0){
    lambda<-pivotVal/A[pivotindex, pivotindex]
  } else
  {
    lambda<-0
  }
  return(lambda)
}

