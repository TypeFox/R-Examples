`BuildOptVec` <-
function(n,m){
  # CAUTION: This 'n' is rather 'nr'
  a1 <- matrix(rep(1,n),ncol=1)
  a2 <- matrix(rep(0,m),ncol=1)
  a3 <- -a1

  rbind(a1,a2,a3)
}

