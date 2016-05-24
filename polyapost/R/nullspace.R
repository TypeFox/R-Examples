#finds the matrix that spans the null space of a known matrix with the no. of columns >= the no. of rows
nullspace<-function(A)
{
  if(! is.matrix(A)) stop("A not matrix")
  n<-nrow(A)
  p<-ncol(A)
  if(n >= p) stop("no. of rows greater or equal to the no. of columns")
  s<-svd(A, nu=n, nv=p)
  k<-0
  for(i in 1:n) {
     if (s$d[i] >1e-6) k<-k+1
  }
  v2<-NULL
  for(i in (k+1):p) v2<-cbind(v2,s$v[,i])
  return(v2)
}
