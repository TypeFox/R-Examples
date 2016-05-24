#' Partitioned Distances
#'
#' Compute distance matrix between two matrices of observations,
#' or two subsets of one matrix
#'
#' @param X a matrix of n observations where columns represent features
#'    of the observations
#' @param Y optional.  A second matrix of p observations like X.  Y
#'    must have the same number of columns as X
#' @param indices.A optional.  A vector of integers representing
#'    row indices from X.  This should only be used when Y is not provided. 
#' @param indices.B optional.  A vector of integers representing
#'    row indices from X.  This should only be used when Y is not provided.
#' @details pdist computes a n by p distance matrix using two seperate
#'    matrices.  pdist allows the user to factor out observations into
#'    seperate matrices to improve computations.  The function dist
#'    computes the distances between all possible pair wise elements,
#'    pdist only computes the distance between observations in X with
#'    observations in Y; distances between observations in X and other
#'    observations in X are not computed, and likewise for Y.
#'
#'    If a second matrix Y is not provided, indices.A and indices.B
#'    can be provided together to specify subsets of X to be computed.
#'    A new matrix X is created by taking X[indices.A,] and Y is
#'    created using X[indices.B,].
#'    
#'    The return value of pdist is a distance vector, much like the default
#'    return value for dist.  However, it can be accessed like a full
#'    distance matrix.  If mypdist = pdist(X,Y), mypdist[i,j] is the
#'    distance between X[i,] and Y[j,].  Similarly, mypdist[i,] is a
#'    vector of distances between X[i,] and all observations in Y.
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.pdist = pdist(x, indices.A = 1:3, indices.B = 8:10)
#'   message("Find the distance between observation 1 and 10 of x")
#'   x.pdist[1,3]
#'   message("Converting a pdist object into a traditional distance matrix")
#'   as.matrix(x.pdist)
#' @export
#' @useDynLib pdist
pdist = function(X, Y = NULL, indices.A = NULL, indices.B = NULL) {
  if (!is.null(Y)) {
    if (identical(X,Y)) warning("Y is the same as X, did you mean to use dist instead?")
  }
  else {
    if (is.null(indices.A) & is.null(indices.B)) {
      warning("Not enough parameters specified: at least Y, or both indices.A and
               indices.B must be specified.  Returning dist(X) instead.")
      return (dist(X))
    }
    if (length(indices.A) == nrow(X) & length(indices.B) == nrow(X)) {
      warning("indices.A and indices.B should be used to make subsets of X.
               This configuration will compute the distance between all pairs
               of rows in X.")
    }
  }
   
  if (is.null(Y) & (!is.null(indices.A) & !is.null(indices.B))) {
      Y = X[indices.B,]
      X = X[indices.A,]
  }
 
  if (is.vector(X)) X = matrix(X, nrow = 1)
  if (is.vector(Y)) Y = matrix(Y, nrow = 1)
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Y)) Y = as.matrix(Y)

  nx = as.integer(nrow(X))
  ny = as.integer(nrow(Y))
  p = as.integer(ncol(X))

  X.vec = as.single(as.vector(X))
  Y.vec = as.single(as.vector(Y))

  distances = single(nrow(X) * nrow(Y))

  result = .C("Rpdist", X.vec, Y.vec, nx, ny, p, distances=distances,
              NAOK = T)
  new("pdist", dist = result$distances, n = nrow(X), p = nrow(Y))
}
