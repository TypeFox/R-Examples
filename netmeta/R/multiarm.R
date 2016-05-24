multiarm <- function(r) {
  ##
  ## Dimension of r and R
  ##
  m <- length(r)                 # Number of edges
  n <- (1 + sqrt(8 * m + 1)) / 2 # Number of vertices
  ##
  ## Construct adjacency matrix and edge.vertex incidence matrix of
  ## complete graph of dimension n
  ##
  A <- 1 - diag(rep(1, n))
  B <- matrix(0, nrow = m, ncol = n)
  k <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      B[k, i] <-  1
      B[k, j] <- -1
    }
  }
  ##
  ## Distribute the edge variances on a symmetrical n x n matrix, R
  ##
  R <- diag(diag(t(B) %*% diag(r) %*% B)) - t(B) %*% diag(r) %*% B
  ##
  ## Construct pseudoinverse Lt from given variance (resistance) matrix R
  ## using a theorem equivalent to Theorem 7 by Gutman & Xiao
  ## Lt <- -0.5 * (R - (R %*% J + J %*% R) / n + J %*%R %*% J / n^2)
  ##
  Lt <- -0.5 * t(B) %*% B %*% R %*% t(B) %*% B / n^2
  ##
  ## Compute Laplacian matrix L from Lt
  ## 
  L <- solve(Lt - 1 / n) + 1 / n
  ##
  ## Compute weight matrix W and variance matrix V from Laplacian L
  ## 
  W <- diag(diag(L)) - L
  V <- 1 / W
  ##
  ## Compute original variance vector v from V
  ##
  v <- rep(0, n)
  edge <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      edge <- edge + 1
      v[edge] <- V[i, j]
    }
  }
  ##
  ## Result
  ##
  res <- list(n = n, r = r, R = R, Lt = Lt, L = L, W = W, V = V, v = v)
  res
}
