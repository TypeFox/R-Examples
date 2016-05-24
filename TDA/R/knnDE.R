knnDE <-
function(X, Grid, k){

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k should be a positive integer")
  }

  d <- ncol(X)
  n <- nrow(X)
  r.k <- apply(FNN::knnx.dist(X, Grid, k = k, algorithm = "kd_tree"), 1, max)
  v.d <- pi^(d/2) /gamma(d/2+1)
  out <- k / (n * v.d * r.k ^ d)  
  return(out)
}