l1pcastar <- function (X, projDim=1, center=TRUE, scores=FALSE, projPoints=FALSE, dispExp=FALSE)
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

  X <- t(X)

  pcLength    <- nrow(X) * nrow(X)
  scoreLength <- 0
  projLength <- 0
  if (scores || dispExp) {
    scoreLength <- projDim * ncol(X)
  }
  if (projPoints) {
    projLength  <- nrow(X) * ncol(X)
  }

  sol <- .C ("l1pcastar", as.double(X), as.integer(dim(X)), as.integer(projDim), as.integer(scores), as.integer(projPoints), loadings=double(pcLength), scores=double(scoreLength), projPoints=double(projLength), PACKAGE="pcaL1")

  solution <- new.env()
  solution$loadings <- matrix(sol[["loadings"]], ncol=nrow(X), byrow=FALSE)
  
  if (scores || dispExp) {
    solution$scores <- matrix(sol[["scores"]], ncol=projDim, byrow=TRUE)
    row.names(solution$scores) <- colnames(X)
    totalDisp <- sum(abs(X))
    scoreDisp <- (apply(abs(solution$scores), 2, sum))
    solution$dispExp <- scoreDisp/totalDisp
  }
  
  if (projPoints) {
    solution$projPoints <- matrix(sol[["projPoints"]], ncol=nrow(X), byrow=TRUE)  
    row.names(solution$projPoints) <- colnames(X)
  }
  
  solution <- as.list(solution)
  class(solution) <- "l1pcastar"
  solution
}
