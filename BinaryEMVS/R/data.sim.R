#' High Dimensional Correlated Data Generation
#'
#' Generates an high dimensional dataset with a subset of columns being related to the response, while
#' controlling the maximum correlation between related and unrelated variables.
#'
#' @param n sample size
#' @param p total number of variables
#' @param pr the number of variables related to the response
#' @param cor the maximum correlation between related and unrelated variables
#'
#' @return Returns an nxp matrix with the first pr columns having maximum correlation cor with
#'         the remaining p-pr columns
#'
#' @examples
#' data=data.sim(n=100,p=1000,pr=10,cor=.6)
#' max(abs(cor(data))[abs(cor(data))<1])
#'
#' @export
data.sim=function(n=100,p=1000,pr=3,cor=.6)
{
  beta.Vec=matrix(0,p,1)
  k=pr    #dimension of diagonal matrix D used to construct T
  #INCLUDING LARGEST EIGENVALUE
  alpha=0  #intercept term=0 w.l.o.g
  
  A = matrix(runif(n*p),nrow=n,ncol=p)
  related=A[,1:pr] 
  
  lamsq=cor^2/(1-cor^2) #largest eigenvalue of T'T
  
  # D is thethe diagonal part of T
  D = diag(c(sqrt(lamsq),sort(runif(k-1,0,sqrt(lamsq)),decreasing=TRUE)))
  
  T=rbind(cbind(D,matrix(0,nrow=k,ncol=n-pr-k)),
          cbind(matrix(0,nrow=pr-k,ncol=k),matrix(0,nrow=pr-k,ncol=n-pr-k)))
  
  #T is block diagonal
  # [ D 0 ]
  # [ 0 0 ]
  
  A.qr <- qr(A)
  Q <- qr.Q(A.qr)
  
  al=Q[,1:pr]     #take the first pr columns of the G-S orthonormalization
  be=Q[,(pr+1):n] #take the remaioning n-pr columns of the G-S orthonormalization
  
  #first construct the linear space that will generate the unrelated p-pr
  #variables, then get the basis for this space
  C=be+al%*%T   
  ro=nrow(C)  #number of columns of this basis, "k" in out notes
  co=ncol(C)
  
  # my explanation for this is the picture I attached in the email.
  b=kronecker(matrix(runif((p-pr)*(co),-.5,.5),p-pr,co),matrix(1,nrow=n,ncol=1))
  C.rep=C[rep(1:n,times=p-pr),]
  variables= apply(b * C.rep,1,sum)
  unrelated=matrix(variables,nrow=n,ncol=p-pr,byrow=F)
  design.Mat = cbind(related, unrelated)
  
  ########## STANDARDIZE DESIGN MATRIX################
  design.Mat=scale(design.Mat,center=TRUE,scale=TRUE)
  
  return(design.Mat)
}