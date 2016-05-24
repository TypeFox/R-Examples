#feasible1 returns a feasible solution for Ax=b, x>=eps, where eps is a small positive number and a negative vector if there is no feasible solution 

feasible1<-function(A, b, eps)
{
  if(! is.matrix(A)) stop("A not a matrix")
  n<-nrow(A)
  p<-ncol(A)
  if(n != length(b)) stop(" no. of rows does not match length of second parameter ")
  if(p != length(eps)) stop("no. of columns  does not match length of third parameter")
  rside<-b-A%*%eps
  s<-0
  for(i in 1:n){
    if(rside[i] >= 0) s<-s+1
  }
  if(s == n){
    d<-diag(rep(1,n))
    linsol<-simplex(a=c(rep(0,p), rep(1,n)), A3=cbind(A,d), b3=rside, maxi=FALSE)
    if(linsol$value<1e-7){
      x0<-as.vector(linsol$soln[1:p])
      x<-x0+eps
    }
    else  x<--1+eps
  }
  else x<--1+eps
  return(x)
}
