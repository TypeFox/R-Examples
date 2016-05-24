library(pcalg)

p <- 10
reps <- 10
set.seed(21)
res <- logical(reps)
for (i in 1:reps) {
  amat <- matrix(sample(c(0,1),p*p,replace=TRUE),p,p)
  diag(amat) <- 0
  skelT <- amat+t(amat)
  skelT[skelT != 0] <- 1

  ## same skeleton?
  my.dag <- amat2dag(amat)
  skelDAG <- my.dag+t(my.dag)
  res1 <- all(skelDAG == skelT)

  ## acyclic?
  res2 <- ggm::isAcyclic(my.dag)

  res[i] <- res1 & res2
}

if(!all(res)) {
  stop("Test amat2dag: Problem")
}
