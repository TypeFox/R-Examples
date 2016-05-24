con.parmsFUN <- function(jk, x, y, n){
  j = jk[1]
  k = jk[2]
  a <- con.parms(x,y,n,j,k,1,1)
  nr <- nrow(a$theta)
  est <- a$theta[nr,  ]
  b<-con.est(x[j],x[k],est)
  s2<-1/b$eta1
  t2<-1/b$eta2
  return(p.ll(n, j, k, s2, t2))
}