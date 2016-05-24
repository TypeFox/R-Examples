getDtfPosSparse <- function(n,ord,pos) {
  D = bandSparse(n, m=n, c(0,1), diagonals=list(rep(-1,n),rep(1,n-1)))
  D0 = D
  for (i in Seq(1,ord)) {
    wts = c(i/(pos[(i+1):n]-pos[1:(n-i)]),rep(1,i))
    D = D0 %*% (wts * D)
  }
  return(D[1:(n-ord-1),])
}
