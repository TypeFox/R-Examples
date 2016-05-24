# Update the QR factorization, after a row has been
# added. Here D1 is (m2+1) x n, Q1 is m2 x min(n,m2),
# and R is min(n,m2) x n.

updateT <- function(y,D1,D2,Q1,Q2,R,q,row) {
  m2 = nrow(D1)-1
  n = ncol(D1)

  # The full D
  D = rbind(D1,D2)
  m1 = nrow(D)
  
  # The full Q, the new Q
  Q = cbind(Q1,Q2)
  Q = rbind(cbind(rep(0,m2),Q),c(1,rep(0,m2)))

  a = .C("update2",
    y=as.double(y),
    D=as.double(D),
    row=as.double(row),
    m=as.integer(m1),
    n=as.integer(n),
    q=as.integer(q),
    dup=FALSE,
    package="genlasso")

  y = a$y
  D = matrix(a$D,nrow=m1)

  # Form D1,Q2 and Q1,Q2 and trim R
  D1 = D[Seq(1,m2+1),,drop=FALSE]
  D2 = D[Seq(m2+2,m1),,drop=FALSE]
  r = min(n,m2+1)
  Q1 = Q[,Seq(1,r),drop=FALSE]
  Q2 = Q[,Seq(r+1,m2+1),drop=FALSE]
  R = rbind(a$row,R)[Seq(1,r),,drop=FALSE]
  
  return(list(y=y,D1=D1,D2=D2,Q1=Q1,Q2=Q2,R=R))
}
