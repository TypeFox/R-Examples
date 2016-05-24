splithalf.r <-
function(x, sims=1000, graph=TRUE, seed=2) {
  N <- ncol(x)
  half.N <- N/2
    if(half.N == round(half.N)) {
      coded <- c(rep(1,half.N),rep(2,half.N))
      }
    if(half.N != round(half.N)) {
      coded <- c(rep(1,half.N+.5),rep(2,half.N-.5))
      }
  
  store <- rep(NA, sims)
  if(seed!=F) {set.seed(seed)}
    for(i in 1:sims) {
      rand.assign <- sample(coded, N, FALSE)
      assign.1 <- x[,rand.assign==1]
      assign.2 <- x[,rand.assign==2]
      comp.1 <- sapply(data.frame(t(assign.1)), meanif, nomiss=0)
      comp.2 <- sapply(data.frame(t(assign.2)), meanif, nomiss=0)
      store[i] <- cor(comp.1, comp.2)
    }

  Avg.r <- fisherz2r(mean(fisherz(store)))
  Up.r <- (2*Avg.r) / (Avg.r+1)
  SD.r <- fisherz2r(sd(fisherz(store)))

  if(graph==T) {
  hist(store, main="Histogram of Split-Half rs", xlab="Split-Half r Values", ylab="Frequency", col="cyan")
  abline(v=Avg.r, col="red")
  }

  out <- cbind(N, Avg.r, Up.r, SD.r)
  colnames(out) = c("N Vars", "Mean Split-Half r", "Rel", "Rel SD")
  rownames(out) = c("Results")
  return(out)
}
