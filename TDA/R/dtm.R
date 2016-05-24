dtm <-
function(X, Grid, m0, weight = 1) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (NCOL(X) != NCOL(Grid)) {
    stop("dimensions of X and Grid do not match")
  }
  if (!is.numeric(m0) || length(m0) != 1 || m0 < 0 || m0 > 1) {
    stop("m0 should be a number between 0 and 1")
  }
  if (!is.numeric(weight) || 
      (length(weight) != 1 && length(weight) != NROW(X))) {
    stop("weight should be either a number or a vector of length equals the number of sample")
  }

  # without weight
  if (length(weight) == 1) {
    X <- as.matrix(X) 
    k0 <- ceiling(m0 * NROW(X))
    distances <- FNN::knnx.dist(X, as.matrix(Grid), k = k0,
        algorithm = c("kd_tree"))
    return (sqrt(apply(distances^2, 1, sum)/k0))

  # with weight
  } else {
    X0 <- as.matrix(X[weight != 0, , drop = FALSE]) 
    weight0 <- weight[weight != 0]
    weight0sort <- sort(weight0)
    weightBound <- m0 * sum(weight0)
    weightSumTemp <- 0
    for (k0 in seq(along = weight0)) {
      weightSumTemp <- weightSumTemp + weight0sort[k0]
      if (weightSumTemp >= weightBound) {
        break
      }
    }
    indexDistance <- FNN::get.knnx(X0, as.matrix(Grid), k = k0,
        algorithm = c("kd_tree"))
    return (Dtm(knnIndex = indexDistance[["nn.index"]],
        knnDistance = indexDistance[["nn.dist"]], weight = weight0,
        weightBound = weightBound))
  }
}
