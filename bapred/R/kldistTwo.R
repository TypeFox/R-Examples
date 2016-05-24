kldistTwo <-
function(xb1, xb2) {

  ##require("FNN")

  distmat1 <- as.matrix(dist(xb1))
  dists1 <- distmat1[lower.tri(distmat1)]

  distmat2 <- as.matrix(dist(xb2))
  dists2 <- distmat2[lower.tri(distmat2)]

  distmat12 <- apply(xb1, 1, function(y) sqrt(rowSums(sweep(xb2, 2, y, "-")^2)))
  dists12 <- as.vector(distmat12)

  (nrow(xb1)*FNN::KL.dist(dists1, dists12, k=5)[5] + nrow(xb2)*
    FNN::KL.dist(dists2, dists12, k=5)[5])/(nrow(xb1) + nrow(xb2))

}
