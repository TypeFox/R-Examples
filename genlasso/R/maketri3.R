# Make the R factor upper triangular by Givens rotating
# its columns. We have to apply the same rotations to y,
# D1,D2. Here D1 is m2 x n, R is m2 x n, and the first
# q columns of D1 and R are zero.

maketri3 <- function(y,D1,D2,R,q) {
  m2 = nrow(D1)
  n = ncol(D1)

  # The full D
  D = rbind(D1,D2)
  m1 = nrow(D)
  
  a = .C("maketri3",
    y=as.double(y),
    D=as.double(D),
    R=as.double(R),
    m1=as.integer(m1),
    m2=as.integer(m2),
    n=as.integer(n),
    q=as.integer(q),
    dup=FALSE,
    package="genlasso")

  y = a$y
  D = matrix(a$D,nrow=m1)
  R = matrix(a$R,nrow=m2)

  # Form D1,D2
  D1 = D[Seq(1,m2),,drop=FALSE]
  D2 = D[Seq(m2+1,m1),,drop=FALSE]
  
  return(list(y=y,D1=D1,D2=D2,R=R))
}
