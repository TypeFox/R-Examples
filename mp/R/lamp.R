#' @importFrom stats dist
{}
#' Local Affine Multidimensional Projection
#'
#' Creates a 2D representation of the data. Requires a subsample
#' (sample.indices) and its 2D representation (Ys).
#'
#' @param X A data frame or matrix.
#' @param sample.indices The indices of data points in X used as subsamples. If
#'   not given, some points from X will be randomly selected and Ys will be generated
#'   by calling forceScheme on them.
#' @param Ys Initial 2D configuration of the data subsamples (will be ignored if
#'   sample.indices is NULL).
#' @param cp Proportion of nearest control points to be used.
#' @return The 2D representation of the data.
#'
#' @references Joia, P.; Paulovich, F.V.; Coimbra, D.; Cuminato, J.A.; Nonato,
#'   L.G., "Local Affine Multidimensional Projection," Visualization and
#'   Computer Graphics, IEEE Transactions on , vol.17, no.12, pp.2563,2571,
#'   Dec. 2011
#'
#' @examples
#' # Iris example
#' emb <- lamp(iris[, 1:4])
#' plot(emb, col=iris$Species)
#'
#' @useDynLib mp
#' @export
lamp <- function(X, sample.indices=NULL, Ys=NULL, cp=1) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  if (is.null(sample.indices)) {
    n <- nrow(X)
    sample.indices <- sample(1:n, sqrt(n))
    Ys <- NULL
  }

  if (is.null(Ys)) {
    sample.indices <- as.vector(sample.indices)
    Ys <- forceScheme(stats::dist(X[sample.indices, ]))
  }

  if (!is.matrix(Ys)) {
    Ys <- as.matrix(Ys)
  }

  if (length(sample.indices) != nrow(Ys)) {
    stop("sample.indices and Ys must have the same number of instances")
  }

  .Call("mp_lamp", X, sample.indices, Ys, cp, PACKAGE="mp")
}

#lamp = function(X, sample.indices = NULL, Ys = NULL) {
#    if (is.null(sample.indices)) {
#        n = nrow(X)
#        sample.indices = sample(1:n, sqrt(n))
#    }
#
#    Xs = X[sample.indices, ]
#
#    if (is.null(Ys)) {
#        Ys = forceScheme(as.matrix(dist(Xs)))
#    }
#
#    Y = t(apply(X, 1, function(point) {
#        alphas = apply(Xs, 1, function(sample.point) 1 / sum((sample.point - point)^2))
#        alphas.sum = sum(alphas)
#        alphas.sqrt = sqrt(alphas)
#
#        X.til = colSums(Xs * alphas) / alphas.sum
#        Y.til = colSums(Ys * alphas) / alphas.sum
#        X.hat = sweep(Xs, 2, X.til, '-')
#        Y.hat = sweep(Ys, 2, Y.til, '-')
#
#        A = X.hat * alphas.sqrt
#        B = Y.hat * alphas.sqrt
#        AtB = t(A) %*% B
#        s = propack.svd(AtB, neig = 2) # We just need the first 2
#        M = s$u %*% s$v
#
#        y = (point - X.til) %*% M + Y.til
#        y
#    }))
#}
