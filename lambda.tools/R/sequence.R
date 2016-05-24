# :vim set filetype=R

#' Safely get an element from a vector
#'
#' This function guarantees a vector of length > 1 as the return value of
#' an indexing operation.
#'
#' @section Usage:
#' item(v, idx)
#'
#' @section Details:
#' Standard R indexing yields different results depending on the input.
#' When either an empty vector or a NULL is passed to the indexing
#' operator, an empty vector is returned. However, if the index is NA,
#' the return value will be a vector of NAs having the same length as
#' the original vector. This inconsistent behavior requires special
#' handling whenever the index value is computed dynamically.
#'
#' This function is designed to create a consistent return value for a 
#' bad index value, which is defined as {NULL, NA, vector of length 0}.
#' If any of these values are used as the index, then NA is returned
#' instead of an empty vector.
#'
#' @name item
#' @param v A sequence
#' @param idx The index of the element to extract
#' @return Either the value of x[idx] or NA for invalid index values
#' 
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' # Compare default behavior with item
#' (1:10)[NA]
#' item(1:10, NA)
#'
#' # Negative indices are still allowed
#' item(1:10, -2) 
item(v, NULL) %as% NA
item(v, NA) %as% NA
item(v, EMPTY) %as% NA
item(v, 0) %as% NA
item(v, idx) %when% { is.null(dim(v)) } %as% v[idx]


#' Pad a vector with some default value
#'
#' This function pads a vector with default values as a way to coerce the
#' value to some predetermined length.
#'
#' @section Usage:
#' pad(x, head, tail=0, default=NA)
#'
#' @section Details:
#' It is common for sequence operations to return a sequence that is
#' shorter than the original sequence. This phenomenon can be
#' annoying when binding the output with the input in a regular
#' data structure like a matrix or data.frame. This function prepends
#' or appends a specified value to a data structure to ensure that the
#' length of the data structure is compatible with another data structure.
#'
#' @name pad
#' @param x A vector to pad
#' @param head The amount to prepend
#' @param tail The amount to append
#' @param default The value to use for the pad
#' @return A padded sequence
#'
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' # A moving average results in n - window + 1 results, so pad at the
#' # head to get a vector of length 50
#' x <- abs(rnorm(50))
#' m <- maprange(x, 10, mean)
#' pad(m, 9) 
#' 
#' # Pad at the end instead of the beginning. Note that the head must
#' # explicitly be set to 0
#' pad(m, 0, 9)
#'
#' # Pad on both sides
#' pad(m, 4, 5) 
#'
#' # Use a different default value
#' pad(m, 9, default=0)
pad(x, head, tail=0, default=NA) %when% {
  is.null(dim(x))
} %as% {
  c(rep(default,head),x, rep(default,tail))
}


#' Remove the head and tail of a data structure
#'
#' Remove the specified number of elements from either the head or
#' tail of a data structure. 
#'
#' @section Usage:
#' chomp(x, head=1, tail=1)
#'
#' @section Details:
#' This function is inspired by the PERL function of the same name. While
#' the PERL version is designed for strings, this version is designed for
#' any indexable data structure, typically containing numbers.
#'
#' @name chomp
#' @param x Any indexable data structure
#' @param head The number of elements to be removed from the head of x
#' @param tail The number of elements to be removed from the tail of x
#' @return A data structure with the head and tail chomped off
#'
#' @author Brian Lee Yung Rowe
#' @seealso \code{\link{pad}}
#'
#' @examples
#' chomp(1:10)
#' chomp(letters)
#'
#' chomp(data.frame(x=1:10, y=1:10), head=2, tail=2)
chomp(x, head=1, tail=1) %when% {
  is.null(dim(x))
  head > 0 
  tail > 0
  head + tail < anylength(x)
} %as% {
  x[(1+head):(length(x)-tail)]
}

chomp(x, head=1, tail=1) %when% {
  head > 0
  tail > 0
  head + tail < anylength(x)
} %as% {
  x[(1+head):(nrow(x)-tail), ]
}


#' Slice a sequence into two adjacent sub-sequences
#'
#' A sequence can be sliced using an explicit pivot point or by using
#' a logical expression.
#'
#' @section Usage:
#' slice(x, pivot, inclusive=FALSE)
#'
#' slice(x, expression)
#'
#' @section Details:
#' This function splits a sequence into two adjacent sub-sequences
#' at a pivot point or based on a logical expression. If a pivot
#' point is chosen, then the inclusive parameter determines whether
#' the value associated with the pivot should be included in both
#' sub-sequences. If FALSE, then the indices of the sub-sequences 
#' will have the form [1, pivot], [pivot + 1, n], where n = |x|. If
#' inclusive is TRUE, then the sub-sequences have indices of
#' [1, pivot], [pivot, n]. Obviously the pivot must be an element
#' of the set of indices of x.
#'
#' An alternative construction is to use an expression to define
#' a slice point. The first sub-sequence corresponds to the
#' values where the expression evaluated to TRUE, while the 
#' second sequence corresponds to values when the expression 
#' evaluated to FALSE.
#'
#' In two dimensions only the first variant of this function is
#' defined, as it cannot be guaranteed that a regular matrix will
#' be generated using an arbitrary expression.
#'
#' @name slice
#' @param x An indexable data structure, typically a vector
#' @param pivot The index of the pivot point in x
#' @param inclusive Whether to include the pivot point in the second 
#'  sub-sequence
#' @param expression A logical expression
#' @return A list containing two sub-sequences or sub-matrices
#'
#' @author Brian Lee Yung Rowe
#'
#' @examples
#' # The number 4 is included in each sub-sequence
#' x <- 1:10
#' slice(x, 4, TRUE)
#'
#' # With expressions, the sub-sequences are not necessarily continguous
#' slice(x, x %% 2 == 0)
#'
#' # Same as above but in two dimensions
#' x <- matrix(1:40, ncol=4)
#' slice(x, 4)
#'
slice(x, pivot, inclusive) %::% a : numeric : logical : list
slice(x, pivot, inclusive=FALSE) %when% {
  is.null(dim(x))
  pivot > 0
  pivot < anylength(x)
} %as% {
  left <- x[1:pivot]
  right <- x[(pivot+as.numeric(!inclusive)):length(x)]
  list(left, right)
}

slice(x, pivot, inclusive=FALSE) %when% {
  pivot < anylength(x)
} %as% {
  left <- x[1:pivot,]
  right <- x[(pivot+as.numeric(!inclusive)):anylength(x),]
  list(left, right)
}


slice(x, expression) %::% a : logical : list
slice(x, expression) %when% {
  is.null(dim(x))
  length(expression) == anylength(x)
} %as% {
  left <- x[expression]
  right <- x[!expression]
  list(left, right)
}

slice(x, expression) %when% {
  length(expression) == anylength(x)
} %as% {
  left <- x[expression,]
  right <- x[!expression,]
  list(left, right)
}


