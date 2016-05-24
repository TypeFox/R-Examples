# Make the R factor upper triangular by Givens rotating
# its rows and columns, appropriately. We have to apply
# the same rotations to y,D1,D2 and Q1,Q2. Here D1 is m2
# x n, and R is m2 x n. Effectively (after we remove the
# zero padding) R is square, and its first q columns are
# zero. The last row of R to have a zero diagonal element
# is k.

maketri4 <- function(y,D1,D2,Q1,Q2,R,q,k) {
  m2 = nrow(D1)
  n = ncol(D1)

  # The full D
  D = rbind(D1,D2)
  m1 = nrow(D)

  # The full Q
  Q = cbind(Q1,Q2)

  # The full R
  if (m2>n) R = rbind(R,matrix(0,m2-n,n))

  a = .C("maketri4",
    y=as.double(y),
    D=as.double(D),
    Q=as.double(Q),
    R=as.double(R),
    m1=as.integer(m1),
    m2=as.integer(m2),
    n=as.integer(n),
    q=as.integer(q),
    k=as.integer(k-1),
    dup=FALSE,
    package="genlasso")

  y = a$y
  D = matrix(a$D,nrow=m1)
  Q = matrix(a$Q,nrow=m2)
  R = matrix(a$R,nrow=m2)

  # Form D1,D2 and Q1,Q2 and trim R
  D1 = D[Seq(1,m2),,drop=FALSE]
  D2 = D[Seq(m2+1,m1),,drop=FALSE]
  r = min(n,m2)
  Q1 = Q[,Seq(1,r),drop=FALSE]
  Q2 = Q[,Seq(r+1,m2),drop=FALSE]
  R = R[Seq(1,r),,drop=FALSE]
  
  return(list(y=y,D1=D1,D2=D2,Q1=Q1,Q2=Q2,R=R))
}
