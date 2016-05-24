l1pca <- function (X, projDim=1, center=TRUE, projPoints=FALSE, initialize="l2pca", tolerance=0.0001, iterations=10)
{
  if (class(X) != "matrix") {
    if (class(X) == "data.frame")
      X <- as.matrix(X)
    else
      X <- matrix(X, ncol = 1)
  }
  if (center) {
    myMedian <- apply(X, 2, median)
    myMedMat <- matrix(rep(myMedian, nrow(X)), ncol = ncol(X), byrow=TRUE)
    X <- X-myMedMat
  }

  if (is.matrix(initialize)) {
    initV <- initialize
  }
  else if (initialize == "l2pca") {
    mypca <- prcomp(X, center=center)
    initV <- mypca$rotation[,1:projDim]
  }
  else if (initialize == "random") {
    initV <- matrix(runif(ncol(X)*projDim), ncol=projDim)
  }
  else {
    # return an error
  }
 

  X <- t(X)

  pcLength    <- projDim * nrow(X)
  scoreLength <- projDim * ncol(X)
 
  
  sol <- .C ("l1pca", as.double(X), as.integer(dim(X)), as.integer(projDim), as.double(tolerance), as.integer(iterations), as.double(initV), loadings=double(pcLength), scores=double(scoreLength), PACKAGE="pcaL1")

  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=projDim, byrow=FALSE)
  
  
  solution$scores <- matrix(sol[["scores"]], ncol=projDim, byrow=FALSE)
  row.names(solution$scores) <- colnames(X)
  totalDisp <- sum(abs(X))
  scoreDisp <- (apply(abs(solution$scores), 2, sum))
  solution$dispExp <- scoreDisp/totalDisp
  
  if (projPoints) {
    solution$projPoints <- as.matrix(t((solution$loadings) %*% t(solution$loadings) %*% X))
  }
  
  mysort            <- sort(solution$dispExp, decreasing=TRUE, index.return=TRUE)
  
  solution$scores <- as.matrix(solution$scores[,mysort$ix])
  solution$loadings <- matrix(solution$loadings[,mysort$ix], ncol=projDim, byrow=FALSE)
  solution$dispExp  <- solution$dispExp[mysort$ix]

  solution <- as.list(solution)
  class(solution) <- "l1pca"
  solution
  
}
