pcal1 <- function (X, projDim=1, center=TRUE, scores=FALSE, projPoints=FALSE, dispExp=FALSE, initialize="l2pca")
{
  if (class (X) != "matrix") {
    if (class (X) == "data.frame")
      X <- as.matrix(X)
    else
      X <- matrix(X, ncol = 1)
  }

  if (center) {
    myMedian <- apply(X, 2, median)
    myMedMat <- matrix(rep(myMedian, nrow(X)), ncol = ncol(X), byrow=TRUE)
    X <- X-myMedMat
  }

  A <- X
  X <- X[apply(abs(X),1,sum) > 0,] # get rid of origin points for algorithm
  X <- t(X)
  pcLength <- projDim * (nrow(X))

  initV <- numeric(ncol(X))
  initMethod <- 0
  if (is.numeric(initialize)) {
    initV <- initialize
  }
  else if (initialize == "maxx") {
    initMethod <- 1
  } 
  else if (initialize == "random") {
    initMethod <- 2
  }
  sol <- .C("pcal1", as.double(X), as.integer(dim(X)), as.integer(projDim), loadings=double(pcLength), as.integer(initMethod), initV, PACKAGE="pcaL1")
  
  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=projDim, byrow=FALSE) 

  if (scores || dispExp) {
    solution$scores <- as.matrix(A %*% solution$loadings)
    totalDisp         <- sum(abs(A))
    scoreDisp         <- (apply(abs(solution$scores),2,sum))
    solution$dispExp <- scoreDisp/totalDisp
  }
  
  if (projPoints) {
    solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
  }

  solution <- as.list(solution)
  class(solution) <- "pcal1"
  solution

}  
