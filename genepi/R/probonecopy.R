probonecopy <-
function(n, hf0, g1idx, g1tbl, hpair) {
  p1copy <- rep(0,n)
  for (i in 1:n) {
    l <- g1idx[i]+1
    nhp <- g1tbl[l,2]
    h1 <- hpair[g1tbl[l,1]+(1:nhp)-1,1] + 1
    h2 <- hpair[g1tbl[l,1]+(1:nhp)-1,2] + 1
    php <- hf0[h1]*hf0[h2]
    p1copy[i] <- php[1]/sum(php)
  }
  p1copy
}

