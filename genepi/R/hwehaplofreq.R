hwehaplofreq <-
function(n, two2k, g1idx, g2idx, g1tbl, hpair, tol=1e-8) {
  hf0 <- hf1 <- rep(0, two2k)
  for (i in 1:n) {
    if (g1idx[i] == 0) {
      hf0[g2idx[i]+1] <- hf0[g2idx[i]+1] + 1
    } else {
      l <- g1idx[i]+1
      nhp <- g1tbl[l,2]
      php <- runif(nhp)
      php <- php/sum(php)
      h1 <- hpair[g1tbl[l,1]+(1:nhp)-1,1] + g2idx[i] + 1
      h2 <- hpair[g1tbl[l,1]+(1:nhp)-1,2] + g2idx[i] + 1
      hf0[h1] <- hf0[h1] + 0.5*php
      hf0[h2] <- hf0[h2] + 0.5*php
    }
  }
  hf0 <- hf0/n
  madhf <- sum(abs(hf0-hf1))
  while (madhf > tol) {
    hf1 <- rep(0, two2k)
    for (i in 1:n) {
      if (g1idx[i] == 0) {
        hf1[g2idx[i]+1] <- hf1[g2idx[i]+1] + 1
      } else {
        l <- g1idx[i]+1
        nhp <- g1tbl[l,2]
        h1 <- hpair[g1tbl[l,1]+(1:nhp)-1,1] + g2idx[i] + 1
        h2 <- hpair[g1tbl[l,1]+(1:nhp)-1,2] + g2idx[i] + 1
        php <- hf0[h1]*hf0[h2]
        php <- php/sum(php)
        hf1[h1] <- hf1[h1] + 0.5*php
        hf1[h2] <- hf1[h2] + 0.5*php
      }
    }
    hf1 <- hf1/n
    madhf <- sum(abs(hf0-hf1))
    hf0 <- hf1
  }
  hf0
}

