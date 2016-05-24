cormax.exch <-
function(n, alpha){
  n.max<- n
  cor.max<- diag(1,n.max)
  cor.max[lower.tri(cor.max)]<- rep(alpha,n.max*(n.max-1)/2)
  cor.max[upper.tri(cor.max)]<- rep(alpha,n.max*(n.max-1)/2)
  return(cor.max)
}
