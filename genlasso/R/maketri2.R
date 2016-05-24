# Make the R factor upper triangular by Givens rotating
# its columns. We have to apply the same rotations to
# y and D. Here D is m x n, R is m x n, and k is the
# order of rank deficiency, i.e. rank(R) = min(m,n)-k.

maketri2 <- function(y,D,R,k) {
  m = nrow(D)
  n = ncol(D)
  
  a = .C("maketri2",
    y=as.double(y),
    D=as.double(D),
    R=as.double(R),
    m=as.integer(m),
    n=as.integer(n),
    k=as.integer(k),
    dup=FALSE,
    package="genlasso")

  y = a$y
  D = matrix(a$D,nrow=m)
  R = matrix(a$R,nrow=m)
  
  return(list(y=y,D=D,R=R))
}
