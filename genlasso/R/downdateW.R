# Downdate the QR factorization, after a column has
# been deleted. Here Q1 is m x n, Q2 is m x k, and
# R is n x n.

downdateW <- function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  
  a = .C("downdate1",
    Q1=as.double(Q1),
    R=as.double(R),
    col=as.integer(col-1),
    m=as.integer(m),
    n=as.integer(n),
    dup=FALSE,
    package="genlasso")

  Q1 = matrix(a$Q1,nrow=m)
  R = matrix(a$R,nrow=n)

  # Re-structure: add a column to Q2, delete one from
  # Q1, and trim R
  Q2 = cbind(Q2,Q1[,n])
  Q1 = Q1[,-n,drop=FALSE]
  R = R[-n,-col,drop=FALSE]

  return(list(Q1=Q1,Q2=Q2,R=R))
}
