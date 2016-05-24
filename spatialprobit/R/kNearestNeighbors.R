# spdep::knearneigh  - k nearest neighbours for spatial weights
# build spatial weight matrix W based on the k nearest neighbors (kNN) (default: k=6)
#
# @param X point coordinates (x, y)
# @param m number of neighbors
# @return sparse matrix (n x n) with nearest neigbors
kNearestNeighbors <- function(x, y, k=6) {
  # number of observations
  n <- length(x)
  
  # spatial weight matrix W based on the k= nearest neighbors
  D <- matrix(NA, n, k)  # (n x k) index matrix to the 6 nearest neigbhors from point i
  for (i in 1:n) {
    px <- x[i]
    py <- y[i]
    # euclidean dist from all points to p
    d <- sqrt((x - px)^2 + (y - py)^2)
    
    # determine the m nearest neighbors (rank 1 is the point itself, ranks 2..(m+1) are the m nearest neighbors
    # TODO: if 2 points have the same distance, e.g. ranks=1, 2, 2.5, 2.5 are odd
    D[i, ] <- which(rank(d) %in% 2:(k+1))
  }
  # sparse matrix representation for spatial weight matrix W
  W <- sparseMatrix(i = rep(1:n, k), j=as.vector(D), x=1/k)
  return(W)
}