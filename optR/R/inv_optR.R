#' Invert a matrix using LU decomposition
#' 
#' @description function invert a matrix A using LU decomposition such that A*inv(A)=I
#' @param A   : Input matrix
#' @param tol : tolerance
#' @return A  : Inverse of Matrix A
#' @export
#' @examples 
#' # Invert the matrix using LU decomposition
#' A<-matrix(c(0.6,-0.4,1, -0.3,0.2,0.5,0.6,-1,0.5), nrow=3,ncol=3, byrow = TRUE) 
#' InvA<-inv.optR(A)
inv.optR<-function(A, tol=1e-7){
  if(!is.matrix(A)) A<-as.matrix(A)
  nCOL<-ncol(A)
  b<-diag(nCOL)
  A<-optR(A, method="LU")
  
  # Solve the inverse Equations
  b<-lapply(1:nCOL, 
        FUN=function(i, A, b){
            y<-forwardsubsitution.optR(A$U, b[A$seq,i])
            c<-optR.backsubsitution(A$U, y)
            return(c)
            }, A, b)
  
  b<-matrix(unlist(b), nrow=nCOL, ncol=nCOL, byrow=FALSE)
  b<-machinePrecision(b) # Machine precision
  return(b)
}


#' Function to address machine precision error
#' 
#' @description function to remove the machine precision error
#' @param A   : Input matrix
#' @return A  : return matrix
machinePrecision<-function(A){
  index<-abs(A)<.Machine$double.eps
  A[index]<-0
  return(A)
}

#' Function to extract Lower and Upper matrix from LU decomposition
#' 
#' @description function to extract Lower and Upper matrix from LU decomposition
#' @param A   : Input matrix
#' @return U  : upper triangular matrix
#' @return L  : Lower triangular matrix
#' @export
#' @examples
#' A<-matrix(c(0,-1,1, -1,2,-1,2,-1,0), nrow=3,ncol=3, byrow = TRUE) 
#' Z<-optR(A, method="LU")
#' LUsplit(Z$U)
LUsplit<-function(A){
  nCOL<-ncol(A)
  U<-matrix(rep(0, each=nCOL*nCOL), nrow = nCOL, byrow=T)
  L<-matrix(rep(0, each=nCOL*nCOL), nrow = nCOL, byrow=T)
  U[upper.tri(A, diag=TRUE)]<-A[upper.tri(A, diag=TRUE)]
  L[lower.tri(A, diag=FALSE)]<-A[lower.tri(A, diag=FALSE)]
  diag(L)<-1
  return(list("U"=U, "L"=L))
}