##***************************************************************************
##' Check designs: size and colnames.
##' 
##' Check that the objects \code{XNew} and \code{X} are compatible
##' for interpolation.
##' 
##' @title Check designs: size and colnames
##' 
##' @docType methods
##' 
##' @aliases checkX,matrix-method checkX,Grid-method 
##'
##' @param X Design or grid object.
##'
##' @param XNew Matrix containing (as rows) the points where
##' interpolation must be done, or numeric vector with length \code{d}
##' where \code{d} is the space dimension provided by \code{X}.
##'
##' @param ... Other arguments to be passed to methods.
##'
##' @return A list with a matrix element \code{XNew} corresponding to
##' the wanted matrix possibly coerced and with colnames matching
##' those of \code{X}.
##'
##' @note The order of the two formal arguments has been changed in
##' version 0.0-4.  Former code using positional matching may lead to
##' unexpected results.
##'
##' @examples
##' X <- matrix(1:12, ncol = 3)
##' colnames(X) <- c("Temp", "Press", "Volume")
##' XNew <- matrix(1:3, ncol = 3)
##' XNewMod <- checkX(X = X, XNew = XNew)$XNew
##' XNewMod
setGeneric("checkX", function(X, XNew, ...) standardGeneric("checkX"))


## if 'XNew' is a vector, transformed it into a matrix
setMethod("checkX",
          signature = c(X = "matrix", XNew = "ANY"),
          definition =
          function(X, XNew) {
              d <- ncol(X)
              if (is.null(dim(XNew))) {
                  if (length(XNew) == d) {
                      XNew <- matrix(XNew, nrow = 1L)
                  } else {
                      stop("'XNew' must be a matrix with ", d, 
                           "columns or a numeric vector of length ", d)
                  }
              } else {
                  XNew <- as.matrix(XNew)
                  if (ncol(XNew) != d) {
                      stop("'XNew' must be a matrix with ", d, 
                           "columns or a numeric vector of length ", d)
                  }
              }
              if (!is.null(nms <- colnames(X))) {
                  if (is.null(colnames(XNew))) {
                      colnames(XNew) <- nms
                  } else if (any(colnames(XNew) != nms)) {
                      stop("'XNew' and 'X' with different colnames")
                  }
              }
              list(X = X, XNew = XNew)    
          })

## Grid case
setMethod("checkX",
          signature = c(X = "Grid", XNew = "ANY"),
          definition = function(X, XNew) {
              d <- dim(X)
              nms <- dimnames(X)
              
              if (is.null(dim(XNew))) {
                  if (length(XNew) == d) {
                      XNew <- matrix(XNew, nrow = 1L)
                      colnames(XNew) <- nms
                  } else {
                      stop("'XNew' must be a matrix with ", d, 
                           "columns or a numeric vector of length ", d)
                  }
              } else {
                  XNew <- as.matrix(XNew)
                  if (ncol(XNew) != d) {
                      stop("'XNew' must be a matrix with ", d, 
                           "columns or a numeric vector of length ", d)
                  }
                  if (is.null(colnames(XNew))) {
                      colnames(XNew) <- nms
                  } else if (any(colnames(XNew) != nms)) {
                      stop("'XNew' and 'X' with different names")
                  }
              }
              
              list(X = X, XNew = XNew)           
          })

##*****************************************************************************
##' Find closest point(s) in a design matrix or object.
##'
##' When \code{X} has class \code{"Grid"} weighting the factors
##' has no impact on the (or a) closest point since for each factor
##' there is normally one (and possibly two) level value(s) which
##' minimises the distance to the corresponding value in the given new
##' point (row) of \code{XNew}.
##' 
##' @title Find closest point(s) in a design or object.
##'
##' @docType methods
##' 
##' @aliases closest,matrix-method closest,Grid-method 
##'
##' @param X A matrix containing design points as rows or a
##' \code{"Grid"} object.
##'
##' @param XNew A vector or matrix containing one or several 'new'
##' design points.
##' 
##' @param ... Other arguments for methods.
##'
##' @return A list with elements
##' \item{index}{
##'
##' Integer vector with length equal to the number of 'new' points. Contains
##' the index of the closest point in \code{X}.
##'
##' }
##' \item{value}{
##'
##' Matrix with one row for each 'new' design point containing the
##' closest point as found in \code{X}.
##'
##' }
##' \item{dist}{
##'
##' Numeric vector with length equal to the number of 'new' points. Contains
##' the minimal distances.
##' 
##' }
##'
##' @examples
##' ## levels at 0, 0.2, 0.4, 0.6, 0.8, 1.0
##' g <- Grid(nlevels = c("x" = 6, "y" = 6, "z" = 6))
##' XNew <- c(0.51, 0.61, 0.23)
##' closest(X = g, XNew = XNew)
##' X <- as.matrix(g)
##' closest(X = X, XNew = XNew)
closest <- function(X, XNew, ...){ }
if (!isGeneric("closest")) {
    setGeneric("closest", function(X, XNew, ...) standardGeneric("closest"))
}

## implementation for X = a matrix 
closestXMat <- function(X, XNew, weights = NULL, ...) {
    XNew <- checkX(X = X, XNew = XNew)$XNew
    ind <- rep(NA, nrow(XNew))
    Dist <- function(x, xnew) mean(abs(x - xnew))
    ds <- rep(NA, nrow(XNew))
    for (i in 1:nrow(XNew)) {
        d <- apply(X, 1, Dist, xnew = XNew[i, ])
        ind[i] <- which.min(d)
        ds[i] <- d[ind[i]]
    }
    list(index = ind,
         value = X[ind, ],
         dist = ds)
}
setMethod("closest",
          signature = c(X = "matrix", XNew = "ANY"),
          definition = closestXMat)

## implementation for X = Grid
closestXGrid <- function(X, XNew, weights = NULL, ...) {
              XNew <- checkX(X = X, XNew = XNew)$XNew
              d <- dim(X)
              nms <- dimnames(X)
              nNew <- nrow(XNew)
              nmsNew <- colnames(XNew)
              if (is.null(weights)) weights <- rep(1, d)
              ds <- matrix(NA, nrow = nNew, ncol = d)
              ind <- matrix(NA, nrow = nNew, ncol = d)
              val <- matrix(NA, nrow = nNew, ncol = d)
              colnames(ind) <- colnames(val) <- colnames(ds) <- dimnames(X)
  
              for (j in 1:d) {
                  xlev <- levels(X)[[j]]
                  out <- outer(XNew[ , j], xlev, function(x, y) abs(x - y))
                  ind[ , j] <- apply(out, 1, which.min)
                  val[ , j] <- xlev[ind[ , j]]
                  ds[ , j] <- apply(out, 1, min)
              }
              
              ds <- ds %*% weights
              
              list(index = X@index[ind],
                   value = val,
                   dist = ds)
          }
setMethod("closest",
          signature = c(X = "Grid", XNew = "ANY"),
          definition = closestXGrid)


##******************************************************************************
##' Round a vector (of levels). 
##'
##' @title Round a vector (of levels)
##'
##' @param x Numeric vector to be rounded.
##'
##' @return Rounded numeric vector.
##' 
##' @note At the time, this function can produce poor results
##' when the relative variation in \code{x} is small.
##'
##' @examples
##' ## good
##' round_levels((1:3) / 10^5)
##' round_levels( 10 + (1:3) / 10^5)
##' ## bad
##' round_levels( 10 + (1:3) / 10^6)
round_levels <- function(x) {
    if (length(x) == 1) return(x)
    mx <- log(min(abs(diff(x))), 10)
    if (mx >= 0) d <- 0
    else d <- -floor(mx)
    round(x, digits = d)
}

##***************************************************************************
##' Range of a \code{Grid} object.
##'
##' @title Range of a \code{Grid} object
##' 
##' @param X An object with S4 class \code{"Grid"}, or an object
##' with S3 class \code{"matrix"} or \code{"data.frame"} that can
##' be coerced into a  \code{"Grid"} object.
##'
##' @return A numeric matrix with two rows representing the min
##' and max of the levels for each dimension of the grid.
##'
##' @examples
##' set.seed(1234)
##' myGD <- randGrid(dim = 4)
##' range_Grid(myGD)
##' range_Grid(as.data.frame(myGD))
##' range_Grid(as.matrix(myGD))
range_Grid <- function(X) {
    if (is.data.frame(X)) {
        if (!all(sapply(X, is.numeric))) {
            stop("all columns muts be numeric")
        }
        X <- as.matrix(X)
    }
    if (is.matrix(X)) {
        res <- apply(X, MARGIN = 2, FUN = range)
        rownames(res) <- c("min", "max")
    } else if (is(X, "Grid")){
        res <- sapply(X@levels, range)
        rownames(res) <- c("min", "max")
    }
    res
}

##***************************************************************************
##' Scale a \code{Grid} object.
##'
##' @title Scale a \code{Grid} object.
##'
##' @param X An object with class \code{"Grid"} or which can be
##' coerced into this class.
##'
##' @param fromRange A numeric vector of length \code{2} (min and max)
##' or a matrix with \code{2} rows (min an max) and one column for
##' each grid dimension in \code{X}. This object gives the original
##' "old" range.
##' 
##' @param toRange A numeric vector of length \code{2} (min and max) or
##' a matrix with \code{2} rows (min an max) and one column for
##' each grid dimension in \code{X}.This object gives the target "new" range
##' replacing the old one.
##' 
##' @return
##' An object with the same class as \code{X} but with the
##' each levels rescaled.
##'
##' @seealso \code{\link{range_Grid}}.
##'
##' @examples
##' myGD <- Grid(levels = list(x = c(1, 3, 10), y = c(0, 1000, 2000, 3000)))
##' scale_Grid(myGD)
##' scale_Grid(as.matrix(myGD))
##' scale_Grid(as.data.frame(myGD))
##' 
scale_Grid <- function(X, fromRange = range_Grid(X), toRange = c(0, 1)) {

    if (is(X, "Grid")) d <- dim(X)
    else if (is(X, "data.frame") || is(X, "matrix")) {
        d <- ncol(X)
    } else {
        stop("'X' must be of S4 class \"Grid\" or be a matrix or a data frame")
    }
    
    if (!is.numeric(fromRange)) stop("'fromRange' must be numeric")
    if (is.matrix(fromRange)) {
        if (nrow(fromRange) != 2L) stop("'fromRange' must have 2 rows")
        fromMin <- fromRange[1L, ]
        fromR <- fromRange[2L, ] - fromRange[1L, ] 
    } else {
        fromMin <- fromRange[1L]
        fromR <- fromRange[2L] - fromRange[1L]
    }
    if (any(fromR <= 0)) stop("invalid 'fromRange'") 
    
    if (!is.numeric(toRange)) stop("'toRange' must be numeric")
    if (is.matrix(toRange)) {
        if (nrow(toRange) != 2L) stop("'toRange' must have 2 rows")
        toMin <- toRange[1L, ]
        toR <- toRange[2L, ] - toRange[1L]
    } else {
        toMin <- toRange[1L]
        toR <- toRange[2L]  - toRange[1L]
    }
    if (any(toR <= 0)) stop("invalid 'toRange'")
    
    fromMin <- rep(fromMin, length.out = d)
    fromR <- rep(fromR, length.out = d)
    toMin <- rep(toMin, length.out = d)
    toR <- rep(toR, length.out = d)

    if (is(X, "Grid")){
        for (i in 1L:d){
            X@levels[[i]] <- toMin[i] +
                toR[i] * (X@levels[[i]] - fromMin[i]) / fromR[i]
        }
    } else if (is(X, "data.frame") || is(X, "matrix")) {
        ## a 'sweep' would be faster...
        for (i in 1L:d){
            X[ , i] <- toMin[i] +
                toR[i] * (X[ , i] - fromMin[i]) / fromR[i]
        }
    }
    return(X)
}

##***************************************************************************
##' Reshape or coerce an array of responses.
##'
##' @title Reshape or coerce an array of responses
##'
##' @param X An object with class \code{"GriData"} or that can be
##' coerced to this class.
##'
##' @param Y A vector of length \code{n} or a matrix with \code{n}
##' rows where \code{n} is the number of nodes in \code{X}. In the
##' first case, \code{Y} is coerced into a matrix with one column.
##'
##' @param dropCol Logical. When \code{Y} is a vector and \code{dropCol}
##' is \code{TRUE}, the dimension of the array is \code{nlevels(X)}, else
##' it is \code{c(nlevels(X), 1)}.
##' 
##' @return An array with dimension \code{c(nlevels(X), ncol(Y))} or
##' \code{nlevels(X)} containing the elements of \code{Y} in
##' correspondence with the levels of \code{X}.
##'
##' @seealso \code{\link{apply_Grid}} to compute a vector of response for
##' a given function.
##'
##' @examples
##' myGD <- Grid(nlevels = c(20, 25))
##' F <- apply_Grid(myGD, branin)
##' aF <- array_Grid(X = myGD, Y = F)
##' 
array_Grid <- function(X, Y, dropCol = TRUE) {

    if (!is(X, "Grid"))  X <- as.Grid(X)
    n <- length(X)
    d <- dim(X)
    
    dimY <- dim(Y)    
    if (!is.null(dimY)) {
        if (dimY[1] != n) {
            stop("when 'Y' is an array, its first dimension",
                 " must be equal to the number of nodes 'n'")
        }
        if (length(dimY) != 2) {
            stop("when 'Y' is an array, it must be a matrix")
        }
        m <- dimY[2]
    }  else { 
        if (length(Y) != n) {
            stop("when 'Y' is a vector, its length",
                 " must be equal to the number of nodes 'n'")
        }
        m <- 1
        dimY <- c(n, 1)
        Y <- array(Y, dim = dimY,
                   dimnames = list(NULL, deparse(substitute(Y))))
    }
    
    L <- lapply(levels(X), round_levels)
    res <- array(NA, dim = c(nlevels(X), m))
    dnm <- list()
    for (i in 1:d) {
        dnm[[i]] <- paste(dimnames(X)[i], L[[i]], sep = "=")
    }
    dimnames(res) <- dnm

    for (j in 1:m) {
        res[slice.index(res, MARGIN = d + 1) == j] <- Y[X@index, j]
    }
    dimnames(res) <- dnm
    if (m == 1 && dropCol) {
        dim(res) <- nlevels(X)
        dimnames(res) <- dnm
    } else {
        dimnames(res) <- c(dnm, paste("resp" = 1L:m)) 
    }
    res

}
