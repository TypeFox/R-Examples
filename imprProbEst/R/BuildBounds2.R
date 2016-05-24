`BuildBounds2` <-
function(n,nr,x.frequency,m){
  b1 <- matrix(x.frequency,ncol=1)/n
  b2 <- 1
  b3 <- -1
  b4 <- matrix(rep(0,nr+m+nr),ncol=1)

  rbind(b1,b2,b3,b4)
}

