##' Search for the enclosing Delaunay convex hull
##' 
##' For \code{t = delaunay(cbind(x, y))}, where \code{(x, y)} is a 2D set of
##' points, \code{tsearch(x, y, t, xi, yi)} finds the index in \code{t}
##' containing the points \code{(xi, yi)}.  For points outside the convex hull
##' the index is \code{NA}.
##' 
##' 
##' @param x X-coordinates of triangluation points
##' @param y Y-coordinates of triangluation points
##' @param t Triangulation, e.g. produced by \code{t = delaunayn(cbind(x, y))}
##' @param xi X-coordinates of points to test
##' @param yi Y-coordinates of points to test
##' @param bary If \code{TRUE} return barycentric coordinates as well as index
##' of triangle.
##' @return If \code{bary} is \code{FALSE}, the index in \code{t} containing
##' the points \code{(xi, yi)}.  For points outside the convex hull the index
##' is \code{NA}. If \code{bary} is \code{TRUE}, a list containing:
##' \item{list("idx")}{the index in \code{t} containing the points \code{(xi,
##' yi)}} \item{list("p")}{a 3-column matrix containing the barycentric
##' coordinates with respect to the enclosing triangle of each point code(xi,
##' yi).}
##' @author David Sterratt
##' @note Based on the Octave function Copyright (C) 2007-2012 David Bateman.
##' @seealso tsearchn, delaunayn
##' @export
tsearch <- function(x, y, t, xi, yi, bary=FALSE) {
  if (!is.vector(x))  {stop(paste(deparse(substitute(x)), "is not a vector"))}
  if (!is.vector(y))  {stop(paste(deparse(substitute(y)), "is not a vector"))}
  if (!is.matrix(t))  {stop(paste(deparse(substitute(t)), "is not a matrix"))}
  if (!is.vector(xi))  {stop(paste(deparse(substitute(xi)), "is not a vector"))}
  if (!is.vector(yi))  {stop(paste(deparse(substitute(yi)), "is not a vector"))}
  if (length(x) != length(y)) {
    stop(paste(deparse(substitute(x)), "is not same length as ", deparse(substitute(y))))
  }
  if (length(xi) != length(yi)) {
    stop(paste(deparse(substitute(xi)), "is not same length as ", deparse(substitute(yi))))
  }
  if (ncol(t) != 3) {
    stop(paste(deparse(substitute(t)), "does not have three columns"))
  }
  storage.mode(t) <- "integer"
  out <- .Call("tsearch", as.double(x), as.double(y), t,
               as.double(xi), as.double(yi), as.logical(bary))
  if (bary) {
    names(out) <- c("idx", "p")
  }
  return(out)
}


##' Search for the enclosing Delaunay convex hull
##' 
##' For \code{t = delaunayn(x)}, where \code{x} is a set of points in \code{d}
##' dimensions, \code{tsearchn(x, t, xi)} finds the index in \code{t}
##' containing the points \code{xi}. For points outside the convex hull,
##' \code{idx} is \code{NA}. \code{tsearchn} also returns the barycentric
##' coordinates \code{p} of the enclosing triangles.
##' 
##' @param x An \code{n}-by-\code{d} matrix.  The rows of \code{x} represent
##' \code{n} points in \code{d}-dimensional space.
##' @param t A \code{m}-by-\code{d+1} matrix. A row of \code{t} contains
##' indices into \code{x} of the vertices of a \code{d}-dimensional simplex.
##' \code{t} is usually the output of delaunayn.
##' @param xi An \code{ni}-by-\code{d} matrix.  The rows of \code{xi} represent
##' \code{n} points in \code{d}-dimensional space whose positions in the mesh
##' are being sought.
##' @param fast If the data is in 2D, use the fast C-based \code{tsearch}
##' function to produce the results.
##' @return A list containing: \item{list("idx")}{An \code{ni}-long vector
##' containing the indicies of the row of \code{t} in which each point in
##' \code{xi} is found.} \item{list("p")}{An \code{ni}-by-\code{d+1} matrix
##' containing the barycentric coordinates with respect to the enclosing
##' simplex of each point in \code{xi}.}
##' @author David Sterratt
##' @note Based on the Octave function Copyright (C) 2007-2012 David Bateman.
##' @seealso tsearch, delaunayn
##' @export
tsearchn <- function(x, t, xi, fast=TRUE) {
  ## Check input
  if (!is.matrix(x))  {stop(paste(deparse(substitute(x)), "is not a matrix"))}
  if (!is.matrix(t))  {stop(paste(deparse(substitute(t)), "is not a matrix"))}
  if (!is.matrix(xi)) {stop(paste(deparse(substitute(xi)), "is not a matrix"))}

  n <- dim(x)[2]                        # Number of dimensions
  if (n==2 && fast) {
    return(tsearch(x[,1], x[,2], t, xi[,1], xi[,2], bary=TRUE))
  }
  nt <- dim(t)[1]                       # Number of simplexes
  m <- dim(x)[1]                        # Number of points in simplex grid
  mi <- dim(xi)[1]                      # Number of points to search for
  ## If there are no points to search for, return an empty index
  ## vector and an empty coordinate matrix
  if (mi==0) {
    return(list(idx=c(), p=matrix(0, 0, n + 1)))
  }
  idx <- rep(NA, mi)
  p <- matrix(NA, mi, n + 1)

  ## Indicies of points that still need to be searched for
  ni <- 1:mi

  degenerate.simplices <- c()
  ## Go through each simplex in turn
  for (i in 1:nt) { 
    ## Only calculate the Barycentric coordinates for points that have not
    ## already been found in a simplex.
    b <- suppressWarnings(cart2bary(x[t[i,],], xi[ni,,drop=FALSE]))
    if (is.null(b)) {
      degenerate.simplices <- c(degenerate.simplices, i)
    } else {

      ## Our points xi are in the current triangle if (all(b >= 0) &&
      ## all (b <= 1)). However as we impose that sum(b,2) == 1 we only
      ## need to test all(b>=0). Note that we need to add a small margin
      ## for rounding errors
      intri <- apply(b >= -1e-12, 1, all)

      ## Set the simplex indicies  of the points that have been found to
      ## this simplex
      idx[ni[intri]] <- i

      ## Set the baryocentric coordinates of the points that have been found
      p[ni[intri],] <- b[intri,]

      ## Remove these points from the search list
      ni <- ni[!intri]

      ## If there are no more points to search for, give up
    if (length(ni) == 0) { break }
    }
  }
  if (length(degenerate.simplices) > 0) {
    warning(paste("Degenerate simplices:", toString(degenerate.simplices)))
  }
  return(list(idx=idx, p=p))
}

##' Conversion of Cartesian to Barycentric coordinates.
##' 
##' Given the Cartesian coordinates of one or more points, compute
##' the barycentric coordinates of these points with respect to a
##' simplex.
##' 
##' Given a reference simplex in \eqn{N} dimensions represented by a
##' \eqn{N+1}-by-\eqn{N} matrix an arbitrary point \eqn{\mathbf{P}} in
##' Cartesian coordinates, represented by a 1-by-\eqn{N} row vector, can be
##' written as \deqn{\mathbf{P} = \mathbf{\beta}\mathbf{X}} where
##' \eqn{\mathbf{\beta}} is a \eqn{N+1} vector of the barycentric coordinates.
##' A criterion on \eqn{\mathbf{\beta}} is that \deqn{\sum_i\beta_i = 1} Now
##' partition the simplex into its first \eqn{N} rows \eqn{\mathbf{X}_N} and
##' its \eqn{N+1}th row \eqn{\mathbf{X}_{N+1}}. Partition the barycentric
##' coordinates into the first \eqn{N} columns \eqn{\mathbf{\beta}_N} and the
##' \eqn{N+1}th column \eqn{\beta_{N+1}}. This allows us to write
##' \deqn{\mathbf{P - X}_{N+1} = \mathbf{\beta}_N\mathbf{X}_N +
##' \mathbf{\beta}_{N+1}\mathbf{X}_{N+1} - \mathbf{X}_{N+1}} which can be
##' written \deqn{\mathbf{P - X}_{N+1} = \mathbf{\beta}_N(\mathbf{X}_N -
##' \mathbf{1}\mathbf{X}_{N+1})} where \eqn{\mathbf{1}} is a \eqn{N}-by-1
##' matrix of ones.  We can then solve for \eqn{\mathbf{\beta}_N}:
##' \deqn{\mathbf{\beta}_N = (\mathbf{P - X}_{N+1})(\mathbf{X}_N -
##' \mathbf{1}\mathbf{X}_{N+1})^{-1}} and compute \deqn{\beta_{N+1} = 1 -
##' \sum_{i=1}^N\beta_i} This can be generalised for multiple values of
##' \eqn{\mathbf{P}}, one per row.
##' 
##' @param X Reference simplex in \eqn{N} dimensions represented by a
##' \eqn{N+1}-by-\eqn{N} matrix
##' @param P \eqn{M}-by-\eqn{N} matrix in which each row is the Cartesian
##' coordinates of a point.
##' @return \eqn{M}-by-\eqn{N+1} matrix in which each row is the
##' barycentric coordinates of corresponding row of \code{P}. If the
##' simplex is degenerate a warning is issued and the function returns
##' \code{NULL}.
##' @author David Sterratt
##' @note Based on the Octave function by David Bateman.
##' @examples
##' ## Define simplex in 2D (i.e. a triangle)
##' X <- rbind(c(0, 0),
##'            c(0, 1),
##'            c(1, 0))
##' ## Cartesian cooridinates of points
##' P <- rbind(c(0.5, 0.5),
##'            c(0.1, 0.8))
##' ## Plot triangle and points
##' trimesh(rbind(1:3), X)
##' text(X[,1], X[,2], 1:3) # Label vertices
##' points(P)
##' cart2bary(X, P)
##' @seealso bary2cart
##' @export
cart2bary <- function(X, P) {
  M <- nrow(P)
  N <- ncol(P)
  if (ncol(X) != N) {
    stop("Simplex X must have same number of columns as point matrix P")
  }
  if (nrow(X) != (N+1)) {
    stop("Simplex X must have N columns and N+1 rows")
  }
  X1 <- X[1:N,] - (matrix(1,N,1) %*% X[N+1,,drop=FALSE])
  if (rcond(X1) < .Machine$double.eps) {
    warning("Degenerate simplex")
    return(NULL)
  }
  Beta <- (P - matrix(X[N+1,], M, N, byrow=TRUE)) %*% solve(X1)
  Beta <- cbind(Beta, 1 - apply(Beta, 1, sum))
  return(Beta)
}

##' Conversion of Barycentric to Cartesian coordinates
##' 
##' Given the baryocentric coordinates of one or more points with
##' respect to a simplex, compute the Cartesian coordinates of these
##' points.
##' 
##' @param X Reference simplex in \eqn{N} dimensions represented by a
##' \eqn{N+1}-by-\eqn{N} matrix
##' @param Beta \eqn{M} points in baryocentric coordinates with
##' respect to the simplex \code{X} represented by a
##' \eqn{M}-by-\eqn{N+1} matrix
##' @return \eqn{M}-by-\eqn{N} matrix in which each row is the
##' Cartesian coordinates of corresponding row of \code{Beta}
##' @examples
##' ## Define simplex in 2D (i.e. a triangle)
##' X <- rbind(c(0, 0),
##'            c(0, 1),
##'            c(1, 0))
##' ## Cartesian cooridinates of points
##' beta <- rbind(c(0, 0.5, 0.5),
##'               c(0.1, 0.8, 0.1))
##' ## Plot triangle and points
##' trimesh(rbind(1:3), X)
##' text(X[,1], X[,2], 1:3) # Label vertices
##' P <- bary2cart(X, beta)
##' points(P)
##' @seealso cart2bary
##' @author David Sterratt
##' @export
bary2cart <- function(X, Beta) {
  return(Beta %*% X)
}
