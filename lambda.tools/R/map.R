# :vim set filetype=R
#' Apply a function over each element of a vector
#'
#' This function implements a map operation over arbitrary indexable
#' data structures. Both 1D and 2D data structures are supported.
#'
#' @section Usage:
#' map(x, fn, acc) \%::\% . : Function : . : .
#'
#' map(x, fn, acc=c())
#'
#' @section Details:
#' While many functions in R are vectorized, some functions only work
#' for scalar input. The map function transforms any scalar-valued
#' function into a vectorized function. This is known as the
#' map-equivalent form of the scalar function.
#'
#' The map operation is implemented for 2D data structures as a
#' column-based operation. If a row-based procedure is desired instead,
#' simply transpose the data structure.
#'
#' Conceptually, the map operation is equivalent to the apply family
#' of functions. The reason for this implementation is primarily for
#' pedagogical purposes.
#'
#' @name map
#' @param x Any indexable data structure
#' @param fn A function applied to each x_i in x
#' @param acc An initial data structure to accumulate the value of f(x_i)
#' @return A sequence representing < f(x_i) > for all x_i in x
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{fold}} \code{\link{maprange}} \code{\link{mapblock}}
#'
#' @references Rowe, Brian Lee Yung. 
#' Modeling Data With Functional Programming In R. Chapman & Hall/CRC Press.
#' Forthcoming.
#'
#' @note This function is implemented using recursion and will 
#' throw an error if the length of \code{x} approaches 
#' \code{getOption('expressions') / 8.0}. This limit
#' is due to R attempting to protect against infinite recursion.
#' See \code{options} for more details.
#'
#' @examples
#' map(-10:10, quantize)
#'
#' # Output a list instead of a vector
#' map(-10:10, quantize, acc=list())
#'
#' # Sum the columns of a matrix
#' map(matrix(1:24, ncol=4), sum)
#'
#' # Sum the columns of a data.frame
#' map(data.frame(a=1:6, b=7:12, c=13:18, d=19:24), sum)
#'
map(x, fn, acc) %::% . : Function : . : .
map(EMPTY, fn, acc) %as% acc

map(x, fn, acc=c()) %when% { 
  is.null(dim(x)) 
} %as% {
  map(x[-1], fn, c(acc, fn(x[[1]])))
}

map(x, fn, acc=c()) %as% {
  map(x[,-1,drop=FALSE], fn, c(acc, fn(x[,1])))
}

#' Apply a function over a rolling range of a data structure
#'
#' Either applies a function over a rolling range of a sequence or
#' multiple sequences bound as a matrix or data.frame.
#'
#' @section Usage:
#' maprange(x, window, fn, do.pad=FALSE)
#' 
#' @section Details:
#' This function is intended to work primarily with time series-like
#' objects where the same statistic is computed over a rolling window
#' of the time series. In other packages this operation is referred to as
#' rollapply (e.g. zoo). This version has two significant differences from 
#' other implementations: 1) it is purely functional, and therefore
#' easy to reason about; 2) it has consistent semantics with the
#' family of map functions.
#'
#' Comparing the code for zoo:::rollapply.zoo, which is close to 100 lines,
#' versus the 3 lines separated into 2 function clauses clearly
#' demonstrates the conciseness inherent in functional programming.
#' Mathematics is known for being very compact and powerful. When
#' utilized appropriately, functional programs share this same property.
#'
#' @name maprange
#' @param x Any indexable data structure
#' @param window The length of the sub-sequence to pass to fn
#' @param fn A function applied to a rolling range of x
#' @param do.pad Whether to pad the output to be the same length as the input
#' @return In the 1D case, a vector of length(x) - window + 1 (unless padded)
#' will be returned. Otherwise a matrix with dimension
#' length(x) - window + 1 by ncol(x) will be returned.
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{map}} \code{\link{mapblock}}
#'
#' @examples
#' # Compute a 5-period moving average over a vector
#' maprange(rnorm(20), 5, mean, do.pad=TRUE)
#'
#' # Same as above, but do it for 4 time series
#' maprange(matrix(rnorm(80),ncol=4), 5, mean, do.pad=TRUE)
maprange(x, window, fn, do.pad=FALSE) %when% {
  is.null(dim(x))
  window <= anylength(x)
} %as% {
  y <- sapply(window:length(x), function(idx) fn(x[(idx-window+1):idx]))
  onlyif(do.pad, function(z) pad(z, window-1), y)
}

maprange(x, window, fn, do.pad=FALSE) %when% {
  ! is.null(dim(x))
  window <= anylength(x)
} %as% {
  sapply(1:ncol(x), function(ydx) maprange(x[,ydx], window, fn, do.pad))
}

#' Apply a function over blocks of a vector
#'
#' This form of map operates on non-overlapping adjacent blocks of
#' a data structure.
#'
#' @section Usage:
#' mapblock(x, window, fn, ...)
#'
#' @section Details:
#' This function is useful primarily in the two-dimensional form. The 
#' use case is when a number of rotation matrices should be applied
#' to a set of points. By collecting all the rotation matrices into
#' a larger matrix, it is easy to produce a map process along the 
#' sub-matrices in a way that doesn't require managing indices.
#'
#' Unlike maprange, mapblock doesn't have a do.pad option. Typical
#' usage scenarios begin by constructing a matrix block that is
#' compatible with some other data structure. Hence given a matrix
#' A with dimensions m x n and a window of length m, it is possible
#' to construct a k x m block matrix B composed of smaller m x m 
#' sub-matrices such that each iteration of mapblock operates on a
#' 1 x m vector against an m x m sub-matrix. The point is that by
#' construction the dimensions must be compatible, so padding
#' after the fact becomes unnecessary.
#'
#' The 1D version is provided for completeness and is equivalent to a
#' 2D map, except on the edge cases.
#'
#' @name mapblock
#' @param x Any indexable data structure
#' @param block The block size used to map over
#' @param fn A function applied to a block
#' @param \dots Optional arguments to pass to sapply
#' @return A vector containing the result of fn applied to each block 
#'
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' # Apply multiple rotation matrices to a set of points
#' a <- matrix(sample(12, 20, replace=TRUE), nrow=2)
#' theta <- 2 * pi * sample(360,4, replace=TRUE) / 360
#' b <- fold(theta, function(d,acc)
#'   cbind(acc,matrix(c(cos(d),sin(d),-sin(d),cos(d)), nrow=2)), c())
#' z <- mapblock(b, 2, function(m) m %*% a, simplify=FALSE)
#'  
#' # The 1D version is equivalent to a 2D map
#' x <- 1:24
#' mapblock(x, 4, sum) == map(matrix(x,nrow=4), sum)
#'
mapblock(x, window, fn, ...) %when% {
  is.null(dim(x))
} %as% {
  s <- seq(1, length(x), by=window)
  sapply(s, function(idx) fn(x[idx:min(length(x), idx+window-1)]), ...)
}

mapblock(x, window, fn, ...) %as% {
  s <- seq(1, ncol(x), by=window)
  sapply(s, function(idx) fn(x[,idx:min(ncol(x), idx+window-1)]), ...)
}

