HWIlr <-function(X,zeroadj=0.5) {
  # computes isometric logratio coordinates for a matrix of genotype counts.
  X[X==0] <- zeroadj
  if(is.matrix(X)) {
     c1 <- (log(X[,1]/X[,3]))/sqrt(2)
     c2 <- (log(X[,1]*X[,3]/(X[,2]^2)))/sqrt(6)
     Y <- cbind(c1,c2)
   }
  if(is.vector(X)) {
    c1 <- (log(X[1]/X[3]))/sqrt(2)
    c2 <- (log(X[1]*X[3]/(X[2]^2)))/sqrt(6)
    Y <- c(c1,c2)
  }
  return(Y)
}
