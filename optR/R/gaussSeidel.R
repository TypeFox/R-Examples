#' Gauss-Seidel based Optimization & estimation
#' 
#' @description Function utilizes the Gauss-Seidel optimization to solve equation Ax=b
#' @param A   : Input matrix
#' @param b   : Response
#' @param x   : Initial solutions
#' @param iter  : Number of Iterations
#' @param tol : Convergence tolerance 
#' @return optimal  : Optimal solutions
#' @return initial  : initial solution
#' @return relaxationFactor : relaxation factor
#' @examples
#' A<-matrix(c(4,-1,1, -1,4,-2,1,-2,4), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(12,-1, 5), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gaussseidel", iter=500, tol=1e-7)
gaussSeidel<-function(A, b, x=NULL, iter=500, tol=1e-7){
  w<-1 # relaxation
  k<-10
  p<-1
  nROW<-nrow(A)
  # Initial solution (random values)
  if(is.null(x)) x<-matrix(rep(0, each=nROW), nrow = nROW, byrow=T)
  xini<-x
  for(i in 1:iter){
    xold<-x
    for(k in 1:nROW){
      x[k]<-(1/A[k,k])*(b[k]-sapply(k, FUN=nonDiagMultipication, A, x))
    }
    dx<-sqrt(t(x-xold)%*%(x-xold))
    if(dx<=tol){
      return(list("optimal"=x, "initial"=xini, "relaxationFactor"=w, "itr.conv"=i))
    } else
    {
      if(i==k){
        dx1<-dx
      } 
      if(i==(k+p)){
        dx2<-dx
        w<-2/(1+sqrt(1-((dx2/dx1)^(1/p))))
      }
    }
  }
  print("Optimization Failed to Converge...")
  return(list("x"=x, "xini"=xini, "relaxationFactor"=w, "itr.conv"=i))
}


#' Non-diagnoal multipication
#' 
#' @description Function for non-diagnoal multipication
#' @param i     : Column Index of Matrix A
#' @param A     : Input matrix
#' @param beta  : Response
#' @return asum: Non-diagnol contribution
nonDiagMultipication<-function(i, A, beta){
  a<-seq(1, nrow(A), by=1)
  asum<-sum(A[a!=i, i]*beta[a!=i])
  return(asum)
}