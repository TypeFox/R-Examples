cv.part <- function(n, k) {
  ntest <- floor(n/k)
  ntrain <- n-ntest
  
  ind <- sample(n)

  trainMat <- matrix(NA, nrow=ntrain, ncol=k)
  testMat <- matrix(NA, nrow=ntest, ncol=k)

  nn <- 1:n
  
  for (j in 1:k) {
    sel <- ((j-1)*ntest+1):(j*ntest)
    testMat[,j] <- ind[sel ]
    sel2 <-nn[ !(nn %in% sel) ]
    trainMat[,j] <- ind[sel2]
  }

  return(list(trainMat=trainMat, testMat=testMat))
}
