`adj.rand.index` <- 
function(group1, group2) {
  a <- length(table(group1))
  N <- length(group1)
  ctab <- matrix(,a, a)
  for (j in 1:a) {
    for (i in 1:a) {
      ctab[j,i] <- length(which(group2[which(group1==i)]==j))
    }
  }
  sumnij <- sum(choose(ctab, 2))
  sumai <- sum(choose(colSums(ctab), 2))
  sumbj <- sum(choose(rowSums(ctab), 2))
  Ntwo <- choose(N, 2)
  ari <- abs((sumnij - (sumai*sumbj)/Ntwo)/(0.5*(sumai+sumbj)-(sumai*sumbj)/Ntwo))
  return(ari)
}

