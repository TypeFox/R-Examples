#' Least-Square Projection
#'
#' Creates a q-dimensional representation of multidimensional data. Requires a
#' subsample (sample.indices) and its qD representation (Ys).
#'
#' @param X A data frame or matrix.
#' @param sample.indices The indices of data points in X used as subsamples. If
#'   not given, some rows from X will be randomly selected and Ys will be generated
#'   by calling forceScheme on them.
#' @param k Number of neighbors used to build the neighborhood graph.
#' @param Ys Initial kD configuration of the data subsamples (will be ignored if
#'   sample.indices is NULL).
#' @param q The target dimensionality.
#' @return The qD representation of the data.
#'
#' @references F. V. Paulovich, L. Nonato, R. Minghim, and H. Levkowitz,
#' Least-Square Projection: A fast high-precision multidimensional projection
#' technique and its application to document mapping, vol. 14, no. 3, pp. 564-575.
#'
#' @examples
#' # Iris example
#' emb <- lsp(iris[, 1:4])
#' plot(emb, col=iris$Species)
#'
#' @useDynLib mp
#' @export
lsp <- function(X, sample.indices=NULL, Ys=NULL, k=15, q=2) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  n <- nrow(X)
  if (is.null(sample.indices)) {
    sample.indices <- sample(1:n, 0.1 * n)
    Ys <- NULL
  }

  if (is.null(Ys)) {
    sample.indices <- as.vector(sample.indices)
    Ys <- forceScheme(dist(X[sample.indices, ]))
  }

  if (!is.matrix(Ys)) {
    Ys <- as.matrix(Ys)
  }

  if (length(sample.indices) != nrow(Ys)) {
    stop("sample.indices and Ys must have the same number of instances")
  }

  if (q != ncol(Ys)) {
    stop("target dimensionality must be the same as Ys'")
  }

  nc <- length(sample.indices)
  A <- matrix(data=0, nrow=n + nc, ncol=n)

  Dx <- as.matrix(dist(X))
  for (i in 1:n) {
    neighbors <- order(Dx[i, ])[1 + 1:k]
    A[i, i] <- 1
    alphas <- Dx[i, neighbors]
    if (any(alphas < 1e-6)) {
      idx <- which(alphas < 1e-6)[1]
      alphas <- 0
      alphas[idx] <- 1
    } else {
      alphas <- 1/alphas
      alphas <- alphas / sum(alphas)
      alphas <- alphas / sum(alphas) # With two we guarantee sum(alphas) == 1
    }
    A[i, neighbors] <- -alphas
  }

  for (i in 1:nc) {
    A[n + i, sample.indices[i]] <- 1
  }

  b <- vector("numeric", n + nc)
  b[1:n] <- 0

  Y <- matrix(data=NA, nrow=n, ncol=q)
  for (j in 1:q) {
    b[n + 1:nc] <- Ys[, j]

    L <- t(A) %*% A
    t <- t(A) %*% b
    S <- chol(L)
    Y[, j] <- backsolve(S, backsolve(S, t, transpose=T))
  }

  Y[sample.indices, ] <- Ys
  Y
}
