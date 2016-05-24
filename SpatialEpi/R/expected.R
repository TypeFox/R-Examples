expected <-
function(population, cases, n.strata){

n <- length(population)/n.strata
E <- rep(0, n)
qNum <- rep(0, n.strata)
qDenom <- rep(0, n.strata)
q <- rep(0, n.strata)


# Compute q: strata-specific rates. Numerator and denominator separately
for(i in 1:n.strata){
  indices <- rep(i, n) + seq(0, n-1)*n.strata
  qNum[i] <- qNum[i] + sum(cases[indices])
  qDenom[i] <- qDenom[i] + sum(population[indices])
}
q <- qNum/qDenom


# Compute E expected counts
for(i in 1:n) {
  indices <- 1:n.strata + (i-1)*n.strata
  E[i] <- sum(population[indices]*q)
}

return(E)
}
