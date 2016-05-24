p.estFUN <- function(jk, x, y, n){
  j = jk[1]
  k = jk[2]
  a <- p.est(x,y,n,j,k)
  s2 <- a$sigma2
  t2 <- a$tau2
  return(p.ll(n, j, k, s2, t2))
}