#' @importFrom stats runif
{}
#' Force Scheme Projection
#'
#' Creates a 2D representation of the data based on a dissimilarity matrix. A few
#' modifications have been made in relation to the method described in the
#' literature: shuffled indices are used to minimize the order dependency
#' factor, only a fraction of delta is used for better stability and a tolerance
#' factor was introduced as a second stop criterion.
#'
#' @param D A dissimilarity structure such as that returned by dist or a full
#'          symmetric matrix containing the dissimilarities.
#' @param Y Initial 2D configuration. A random configuration will be used when
#'          omitted.
#' @param max.iter Maximum number of iterations that the algorithm will run.
#' @param tol The tolerance for the accumulated error between iterations. If set
#'            to 0, the algorithm will run max.iter times.
#' @param fraction Controls the point movement. Larger values means less
#'                 freedom to move.
#' @param eps Minimum distance between two points.
#' @return The 2D representation of the data.
#'
#' @references Eduardo Tejada, Rosane Minghim, Luis Gustavo Nonato: On improved
#'   projection techniques to support visual exploration of multi-dimensional
#'   data sets. Information Visualization 2(4): 218-231 (2003)
#'
#' @examples
#' # Eurodist example
#' emb <- forceScheme(eurodist)
#' plot(emb, type = "n", xlab ="", ylab ="", asp=1, axes=FALSE, main="")
#' text(emb, labels(eurodist), cex = 0.6)
#'
#' # Iris example
#' emb <- forceScheme(dist(iris[,1:4]))
#' plot(emb, col=iris$Species)
#' @seealso \code{\link[stats]{dist}} (stats) and \code{\link[proxy]{dist}}
#'   (proxy) for d computation
#'
#' @useDynLib mp
#' @export
forceScheme <- function(D, Y=NULL, max.iter=50, tol=0, fraction = 8.0, eps = 1e-5) {
  # convert d to a matrix
  D <- as.matrix(D)
  if (!is.symmetric(D)) {
    stop("distances matrix is not symmetric")
  }

  n <- nrow(D)
  if (is.null(Y)) {
    Y <- matrix(stats::runif(n*2, 0, 1), ncol=2)
  }

  # check initial configuration
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (nrow(Y) != n) {
    stop("initial configuration does not match dataset size");
  }

  .Call("mp_forceScheme", D, Y, max.iter, tol, fraction, eps, PACKAGE = "mp")
}

