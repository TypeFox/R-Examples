# Update the QR factorization, after a column has been
# added. Here Q1 is m x n, Q2 is m x k, and R is n x n.

updateW <- function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  k = ncol(Q2)
  
  a = .C("update1",
    Q2=as.double(Q2),
    w=as.double(t(Q2)%*%col),
    m=as.integer(m),
    k=as.integer(k),
    dup=FALSE,
    package="genlasso")

  Q2 = matrix(a$Q2,nrow=m)
  w = c(t(Q1)%*%col,a$w)

  # Re-structure: delete a column from Q2, add one to
  # Q1, and expand R
  Q1 = cbind(Q1,Q2[,1])
  Q2 = Q2[,-1,drop=FALSE]
  R = rbind(R,rep(0,n))
  R = cbind(R,w[Seq(1,n+1)])

  return(list(Q1=Q1,Q2=Q2,R=R))
}
