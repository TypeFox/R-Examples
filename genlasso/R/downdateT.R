# Downdate the QR factorization, after a row has been
# deleted. Here Q1 is m x min(n,m) and R is min(n,m) x n.

downdateT <- function(Q1,Q2,R,row) {
  m = nrow(Q1)
  n = ncol(R)
  
  # The full Q
  Q = cbind(Q1,Q2)
  
  # Rearrange Q so that the first row is the one that
  # should to be removed
  Q = rbind(Q[row,],Q[-row,])

  # The full R
  if (m>n) R = rbind(R,matrix(0,m-n,n))

  a = .C("downdate2",
    Q=as.double(Q),
    R=as.double(R),
    m=as.integer(m),
    n=as.integer(n),
    dup=FALSE,
    package="genlasso")

  Q = matrix(a$Q,nrow=m)
  R = matrix(a$R,nrow=m)

  # Form Q1,Q2 and trim R
  r = min(m,n+1)
  Q1 = Q[Seq(2,m),Seq(2,r),drop=FALSE]
  Q2 = Q[Seq(2,m),Seq(r+1,m),drop=FALSE]
  R = R[Seq(2,r),,drop=FALSE]
  
  return(list(Q1=Q1,Q2=Q2,R=R))
}
