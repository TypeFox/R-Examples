#description of the function probvect
#It generates as many points as we want from the polytope A1x=b1, A2x<=b2. It returns #every step-th point.
#length of chain=n*step
probvect<-function(A1, A2, nrow, b2, initsol, step, n)
{
  if(! is.matrix(A1)) stop("A1 not a matrix")
  if(! is.matrix(A2)) stop("A2 not a matrix")
  if(nrow(A2)!=nrow) {stop("wrong no. of rows")}
  P<-nullspace(A1)%*%t(nullspace(A1))
  mat<-NULL
  for(i in 1:n){
    initsol<-probvect1(P, ncol(P), A2, nrow, b2, initsol, step)
    mat<-rbind(mat, initsol)
  }
  return(mat)
}
