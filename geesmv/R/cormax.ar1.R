cormax.ar1 <-
function(n, alpha){
  n.max<- max(n)
  cor.max<- diag(1,n.max)
  lowertri<- rep(0,0)
  for(j in (n.max-1):1){
    lowertri<- c(lowertri,1:j)
  }
  cor.max[lower.tri(cor.max)]<- alpha^lowertri
  cor.max[upper.tri(cor.max)]<- alpha^lowertri[length(lowertri):1]
  return(cor.max)
}
