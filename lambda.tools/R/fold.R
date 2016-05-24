# :vim set filetype=R
#' Successively apply a function to the elements of a sequence 
#'
#' Apply a function to each element of a sequence and the accumulated value
#' of the previous function applications
#'
#' @section Usage:
#' fold(x, fn, acc, ...) \%::\% . : Function : . : ... : .
#'
#' fold(x, fn, acc, ...)
#'
#' @section Details:
#' The fold operation is a generalization of the summation and product
#' operators in mathematics. The idea is that the elements of a sequence
#' can have a function applied to them and then can be aggregated in
#' some arbitrary way. In terms of the summation operator, the general
#' structure is sum f(x_i). This means that the function f is applied
#' to each element of x and then added to some intermediate accumulator.
#' This is equivalent to a function f' : A x B -> B where the single
#' function is responsible for both applying f and also aggregating
#' the accumulated value.
#'
#' A 2D fold is similar to a 2D map in the sense that the function
#' operates on the columns of x. This indicates that fn takes a
#' vector and not a scalar as the first argument. If fn is vectorized,
#' then the behavior of fold will be equivalent to a 2D map over the rows!
#'
#' @name fold
#' @param x Any indexable data structure
#' @param fn A binary operator
#' @param acc Accumulator
#' @return An object containing the accumulated result.
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{map}} \code{\link{foldrange}} \code{\link{foldblock}}
#' @references Haskell Wiki, http://www.haskell.org/haskellwiki/Fold
#' @references Brian Lee Yung Rowe, 
#' Modeling Data With Functional Programming In R.
#'
#' @examples
#' x <- 1:10
#'
#' # This is equivalent to the summation operator
#' sum(x) == fold(x, function(a,b) a+b, 0)
#' sum(x^2) == fold(x, function(a,b) a^2 + b, 0)
#'
#' # This is equivalent to the product operator
#' prod(x) == fold(x, function(a,b) a*b, 1)
#'
#' # Note the equivalence with map
#' x <- matrix(1:24, ncol=4)
#' map(t(x), function(a) sum(a)) == fold(x, function(a,b) a + b, 0)
#'
fold(x, fn, acc) %::% . : Function : . : .
fold(EMPTY, fn, acc) %as% acc

fold(x, fn, acc, ...) %::% . : Function : . : ... : .
fold(x, fn, acc, ...) %when% { 
  is.null(dim(x))
} %as% {
  #fold(x[-1], fn, fn(x[[1]], acc))
  sapply(x, function(xi) {
    acc <<- fn(xi, acc)
    NULL
  }, ...)
  acc
}

fold(x, fn, acc, ...) %as% { 
  #fold(x[,-1,drop=FALSE], fn, fn(x[,1], acc))
  apply(x, 2, function(xi) {
    acc <<- fn(xi, acc)
    NULL
  }, ...)
  acc
}

#' Successively apply a function to a rolling range of a sequence
#'
#' Apply a function to a rolling range of a sequence and the
#' accumulated value of the previous function applications
#'
#' @section Usage:
#' foldrange(x, window, fn, acc, idx) %::% . : numeric : Function : . : numeric : .
#' foldrange(x, window, fn, acc, 0)
#' foldrange(x, window, fn, acc=0, idx=length(x)-window+1)
#'
#' When ! is.null(dim(x))
#' foldrange(x, window, fn, acc=0)
#'
#' @section Details:
#' This function is the fold counterpart of maprange. It's primarily
#' here for completeness purposes, as the utility of this function
#' is still to be determined.
#'
#' @name foldrange
#' @param x Any indexable data structure
#' @param window The length of the sub-sequence passed to fn
#' @param fn The function applied to the rolling range 
#' @param acc An object that stores the intermediate accumulated result
#' @return The accumulated result
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{map}} \code{\link{fold}} \code{\link{foldblock}}
#'
#' @examples
#' \dontrun{
#' # Mean of rolling means
#' z <- sapply(1:500, 
#'   function(x) foldrange(rnorm(50), 10, function(a,b) mean(a) + b) / 41)
#' }
#'
foldrange(x, window, fn, acc, idx) %::% . : numeric : Function : . : numeric : .
foldrange(x, window, fn, acc, 0) %as% acc

foldrange(x, window, fn, acc=0, idx=length(x)-window+1) %when% {
  is.null(dim(x))
  window < anylength(x)
} %as% {
  foldrange(x, window, fn, fn(x[idx:(idx+window-1)], acc), idx-1)
}

foldrange(x, window, fn, acc=0) %when% {
  ! is.null(dim(x))
  window < anylength(x)
} %as% {
  sapply(1:ncol(x), function(ydx) foldrange(x[,ydx], window, fn, acc))
}

#' Successively apply a function to adjacent blocks of a sequence
#'
#' Apply a function to non-overlapping sub-sequences and the 
#' accumulated value of the function application 
#'
#' @section Usage:
#' foldblock(x, window, fn, acc=0)
#'
#' @section Details:
#' This function is the fold counterpart of mapblock. Like mapblock
#' the usefulness of this function is for the 2D case, as it can
#' simplify interacting with matrices. See the example below for
#' using foldblock as a summation operator over matrices
#'
#' @name foldblock
#' @param x Any indexable data structure
#' @param window The number of elements in each sub-sequence
#' @param fn The function applied to the sub-sequence
#' @param acc The intermediate accumulated value
#' @return The accumulated value
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{map}} \code{\link{fold}} \code{\link{foldrange}}
#'
#' @examples
#' # Sum 5 2 x 2 matrices
#' ms <- matrix(sample(40,20, replace=TRUE), nrow=2)
#' foldblock(ms,2, function(a,b) a + b)
#'
#' # 1D foldblock is equivalent to 2D fold
#' x <- 1:12
#' f <- function(a,b) mean(a) + b
#' foldblock(x,3,f) == fold(matrix(x, nrow=3),f, 0)
#'
foldblock(x, window, fn, acc) %::% . : numeric : Function : . : .
foldblock(EMPTY, window, fn, acc) %as% acc

foldblock(x, window, fn, acc=0) %when% { 
  is.null(dim(x))
} %as% {
  y <- min(length(x),window)
  foldblock(x[-(1:y)], window, fn, fn(x[1:y], acc))
}

foldblock(x, window, fn, acc=0) %as% {
  y <- min(ncol(x),window)
  foldblock(x[,-(1:y)], window, fn, fn(x[,1:y], acc))
}
