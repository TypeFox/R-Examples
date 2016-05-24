sum.intdiff.lgamma <-
function(G1,G2){
  l <- length(G1)
  RES <- sapply(1:l,function(x) if (G2[x] - G1[x] > 0) {sum(log(G1[x]:(G2[x]-1)))} else 0)
  return(sum(RES))
}
