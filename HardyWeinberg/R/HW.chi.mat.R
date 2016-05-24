HW.chi.mat <- function(X) {
  n <- rowSums(X)
  nA <- 2*X[,1] + X[,2]
  nB <- 2*X[,3] + X[,2]
  chi <- (4*X[,1]*X[,3] - X[,2]*X[,2])/(nA*nB)
  chi <- n*chi*chi
}

