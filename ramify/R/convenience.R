#' Row/Column Max/Min Indices
#'
#' Returns the indices of the maximum or minimum values along an axis.
#' 
#' @param x A matrix.
#' @param rows If \code{TRUE} (the default) the indices of each row max/min is
#'             returned. 
#' @return A vector of indices.
#' @export
#' @examples
#' m <- mat("94, 20, 44; 40, 92, 51; 27, 69, 74")
#' argmax(m)
#' argmin(m)
argmax <- function(x, rows = TRUE) {
  if (rows) {
    apply(x, MARGIN = 1, which.max)
  } else {
    apply(x, MARGIN = 2, which.max)
  }
}

#' @rdname argmax
#' @export
argmin <- function(x, rows = TRUE) {
  if (rows) {
    apply(x, MARGIN = 1, which.min)
  } else {
    apply(x, MARGIN = 2, which.min)
  }
}


#' View Input as an Array with at Least Two Dimensions.
#' 
#' Ensure that the input has at least two dimensions.
#' 
#' @param x An R object, for example a vector, matrix, array, or data frame.
#' @return The same object, but with a \code{"dim"} attribute.
#' @export
#' @examples
#' x <- 1:10
#' x
#' atleast_2d(x)
atleast_2d <- function(x) {
  if (is.null(dim(x))) {
    dim(x) <- c(length(x), 1L)
  }
  x
}

#' Identity Matrix
#' 
#' Creates an \code{nrow}-by-\code{ncol} identity matrix.
#' 
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @return A \code{nrow}-by-\code{ncol} identity matrix of class 
#'         \code{c("matrix", "mat")}.
#' @seealso \code{\link{diag}}.
#' @export
#' @examples
#' eye(4)  # 4-by-4 identity matrix
#' eye(4, 4)  # 4-by-4 identity matrix
#' eye(3, 5)  # 3-by-5 identity matrix
#' eye(5, 3)  # 5-by-3 identity matrix
eye <- function(nrow = 1, ncol = nrow) {
  diag(1L, nrow, ncol)
}


#' Fill a Matrix
#'
#' Create a matrix filled with the value \code{x}.
#' 
#' @param x The (single) value to fill the matrix with.
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @param ... Further dimensions of the array.
#' @param atleast_2d Logical indicating whether or not to force column vectors 
#'   to have a second dimension equal to one. Defaults to \code{FALSE}. This 
#'   behavior can also be changed globally using, for example 
#'   \code{options(atleast_2d = TRUE)}.
#' @return A matrix or array filled with the value \code{x}. If the result is a
#'         matrix, it will have class \code{c("matrix", "mat")}.
#' @seealso \code{\link{ones}}, \code{\link{zeros}}, \code{\link{falses}}, 
#'   \code{\link{trues}}, \code{\link{mat}}, \code{\link{matrix}}.
#' @export
#' @examples
#' fill(pi, 3, 5)  # 3-by-5 matrix filled with the value of pi
#' fill(pi, 3, 5, 2, 2)  # 3-by-5-by-2-by-2 array filled with the value of pi
#' pi * ones(3, 5)
#' zeros(10)
#' zeros(10, atleast_2d = TRUE)
fill <- function(x, nrow = 1, ncol = 1, ..., atleast_2d = NULL) {
  if (length(list(...)) == 0) {
    if (is.null(atleast_2d)) atleast_2d <- getOption("atleast_2d")
    if (!atleast_2d && ncol == 1) {  # no dim attribute!
      rep_len(x, length.out = nrow)
    } else {
      matrix(x, nrow = nrow, ncol = ncol)
    }
  } else {
    array(x, dim = c(nrow, ncol, unlist(list(...))))
  }
}


#' @rdname fill
#' @export
falses <- function(nrow = 1, ncol = 1, ..., atleast_2d = NULL) {
  fill(FALSE, nrow = nrow, ncol = ncol, ..., atleast_2d = atleast_2d)
}

#' @rdname fill
#' @export
trues <- function(nrow = 1, ncol = 1, ..., atleast_2d = NULL) {
  fill(TRUE, nrow = nrow, ncol = ncol, ..., atleast_2d = atleast_2d)
}

#' @rdname fill
#' @export
ones <- function(nrow = 1, ncol = 1, ..., atleast_2d = NULL) {
  fill(1L, nrow = nrow, ncol = ncol, ..., atleast_2d = atleast_2d)
}

#' @rdname fill
#' @export
zeros <- function(nrow = 1, ncol = 1, ..., atleast_2d = NULL) {
  fill(0L, nrow = nrow, ncol = ncol, ..., atleast_2d = atleast_2d)
}


#' Flatten Matrices/Arrays
#'
#' Flatten (i.e., collapse) a matrix or array to one dimension.
#' 
#' @param x A matrix object.
#' @param across Character string specifying whether to flatten the matrix 
#'   across \code{"rows"} (default) or \code{"columns"}. This option is ignored
#'   for multi-way arrays.
#' @return A numeric vector.
#' @seealso \code{\link{mat}}.
#' @export
#' @examples
#' m <- mat("2, 4, 6, 8; 10, 12, 14, 16")
#' flatten(m)
#' flatten(m, across = "columns")
flatten <- function(x, across = c("rows", "columns")) {
  if (is.matrix(x)) {
    across <- match.arg(across)
    if (across == "rows") x <- t(x)
  } 
  dim(x) <- NULL  # remove dimension attribute
  x
}


#' Matrix Inverse
#' 
#' Calculates the inverse of a square matrix.
#' 
#' @param x A square numeric or complex matrix
#' @param ... Additional optional arguments.
#' @seealso \code{\link{solve}}.
#' @details See the documentation for the \code{base} function 
#'         \code{\link{solve}}.
#' @export
#' @examples
#' m <- 3 * eye(5)
#' inv(m)
inv <- function(x, ...) {
  if (!is.matrix(x)) {
    stop('Argument should be a matrix.', call. = FALSE)
  }
  if (dim(x)[1L] != dim(x)[2L]) {
    stop('Argument should be a square matrix.', call. = FALSE)
  }
  b <- diag(1, nrow(x))
  colnames(b) <- rownames(x)
  solve(x, b, ...)

}


#' Concatenate Matrices
#' 
#' Concatenate matrices along the first or second dimension.
#' 
#' @param ... Vectors or matrices.
#' @return A matrix of class \code{c("matrix", "mat")} formed by combining the
#'          \code{...} arguments column-wise (\code{hcat}) or row-wise 
#'          (\code{vcat}).
#' @seealso \code{\link{bmat}}, \code{\link{cbind}}, \code{\link{rbind}}.
#' @export
#' @examples
#' m1 <- mat("1, 2, 3; 4, 5, 6")
#' m2 <- mat("7, 8, 9; 10, 11, 12")
#' hcat(m1, m2)  # same as 'bmat("m1, m2")'
#' vcat(m1, m2)  # same as 'bmat("m1; m2")'
hcat <- function(...) {
  do.call(cbind, list(...))
}

#' @rdname hcat
#' @export
vcat <- function(...) {
  do.call(rbind, list(...))
}


#' Linearly-spaced Elements
#' 
#' Construct a vector of \code{n} linearly-spaced elements from \code{a} 
#' to \code{b}. 
#' 
#' @param a The starting value of the sequence.
#' @param b The final value of the sequence.
#' @param n The number of samples to generate. Default is 50.
#' @return A vector of linearly-spaced elements.
#' @seealso \code{\link{logspace}}, \code{\link{seq}}.
#' @export
#' @examples
#' linspace(0, 1)
#' linspace(1, 5, 5)
#' linspace(1+2i, 10+10i, 8)
#' logspace(0, pi, 10)
linspace <- function(a, b, n = 50) {
  seq(from = a, to = b, length.out = n)
}


#' Logarithmically-spaced Elements
#' 
#' Construct a vector of \code{n} logarithmically-spaced elements from 
#' 10^\code{a} to 10^\code{b}. 
#' 
#' @param a \code{base^a} is the starting value of the sequence.
#' @param b \code{base^b} is the final value of the sequence.
#' @param n The number of samples to generate. Default is 50.
#' @param base The base of the log space.
#' @return A vector of logarithmically-spaced elements.
#' @note
#' If \code{b = pi} and \code{base = 10}, the points are between 
#' \code{10^a} and \code{pi}, not \code{10^a} and \code{10^pi}, for 
#' compatibility with the corresponding MATLAB/Octave, and NumPy functions.
#' @seealso \code{\link{linspace}}, \code{\link{seq}}.
#' @export
logspace <- function(a, b, n = 50, base = 10) {
  if (b == pi && base == 10)  {
    b <- log(b, base = base)
  }
  base ^ seq(from = a, to = b, length.out = n)
}


#' Rectangular 2-D Grid
#' 
#' Creates matrices for vectorized evaluations of 2-D scalar/vector fields over 
#' 2-D grids.
#' 
#' @param x Numeric vector representing the first coordinate of the grid.
#' @param y Numeric vector representing the second coordinate of the grid.
#' @return a list of matrices, each having class \code{c("matrix", "mat")}.
#' @seealso \code{\link{expand.grid}}, \code{\link{outer}}.
#' @export
#' @examples
#' mg <- meshgrid(linspace(-4*pi, 4*pi, 27))  # list of input matrices
#' z <- cos(mg[[1]]^2 + mg[[2]]^2) * exp(-sqrt(mg[[1]]^2 + mg[[2]]^2)/6)
#' image(z, axes = FALSE)  # color image
#' contour(z, add = TRUE, drawlabels = FALSE)  # add contour lines
meshgrid <- function(x, y = x) {
  lenx <- length(x)
  leny <- length(y)
  list(matrix(rep(x, each = leny), nrow = leny, ncol = lenx),
       matrix(rep(y, times = lenx), nrow = leny, ncol = lenx))
  #   list(as.mat(matrix(rep(x, each = leny), nrow = leny, ncol = lenx)),
  #        as.mat(matrix(rep(y, times = lenx), nrow = leny, ncol = lenx)))
}


#' Matrix/Array of Uniform Random Numbers
#' 
#' Construct a matrix or multi-way array of uniform random deviates.
#' 
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @param ... Further dimensions of the array.
#' @param min Lower limit for the uniform distribution. Must be finite. 
#'   (\code{rand} only).
#' @param max Upper limit for the uniform distribution. Must be finite. 
#'   (\code{rand} only).
#' @param atleast_2d Logical indicating whether or not to force column vectors 
#'   to have a second dimension equal to one. Defaults to \code{FALSE}. This 
#'   behavior can also be changed globally using, for example 
#'   \code{options(atleast_2d = TRUE)}.
#' @return A  matrix or array of pseudorandom numbers. If the result is a
#'         matrix, it will have class \code{c("matrix", "mat")}.
#' @seealso \code{\link{randi}}, \code{\link{randn}}, \code{\link{runif}}.
#' @export
#' @importFrom stats runif
#' @examples
#' rand(100, 100)  # 100 by 100 matrix of uniform random numbers
#' rand(2, 3, min = 100, max = 200)  
rand <- function(nrow = 1, ncol = 1, ..., min = 0, max = 1, atleast_2d = NULL) {
  if (length(list(...)) == 0) {
    if (is.null(atleast_2d)) atleast_2d <- getOption("atleast_2d")
    if (!atleast_2d && ncol == 1) {  # no dim attribute!
      runif(nrow, min = min, max = max)
    } else {
      matrix(runif(nrow * ncol, min = min, max = max), nrow = nrow, 
             ncol = ncol)
    }
  } else {
    array(runif(nrow * ncol * prod(unlist(list(...))), min = min, max = max), 
          dim = c(nrow, ncol, unlist(list(...))))
  }
  
}


#' Matrix/Array of Uniform Random Integers
#' 
#' Construct a matrix or multi-way array of uniform random integers.
#' 
#' @param imax A positive integer.
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @param ... Further dimensions of the array.
#' @param atleast_2d Logical indicating whether or not to force column vectors 
#'   to have a second dimension equal to one. Defaults to \code{FALSE}. This 
#'   behavior can also be changed globally using, for example 
#'   \code{options(atleast_2d = TRUE)}.
#' @return A matrix or array of pseudorandom numbers. If the result is a matrix, 
#'         it will have class \code{c("matrix", "mat")}.
#' @seealso \code{\link{rand}}, \code{\link{randn}}, \code{\link{sample}}.
#' @export
#' @examples
#' randi(2, 5, 5)
randi <- function(imax, nrow, ncol = 1, ..., atleast_2d = NULL) {
  if (!is.integer(imax)) imax <- as.integer(imax)
  if (imax < 1)  {  # make sure imax is a positive integer
    stop("imax must be a positive integer.") 
  }
  if (length(list(...)) == 0) {
    if (is.null(atleast_2d)) atleast_2d <- getOption("atleast_2d")
    if (!atleast_2d && ncol == 1) {  # no dim attribute!
      sample(imax, size = nrow, replace = TRUE)
    } else {
      matrix(sample(imax, size = nrow * ncol, replace = TRUE), nrow = nrow, 
             ncol = ncol)
    }
  } else {
    array(sample(imax, size = nrow * ncol * prod(unlist(list(...))), 
                 replace = TRUE), dim = c(nrow, ncol, unlist(list(...))))
  }
}


#' Matrix/Array of Normal Random Numbers
#' 
#' Construct a matrix or multi-way array of normal random deviates.
#' 
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @param ... Further dimensions of the array.
#' @param mean Mean for the normal distribution. (\code{randn} only).
#' @param sd Standard deviation for the normal distribution. 
#'   (\code{randn} only).
#' @param atleast_2d Logical indicating whether or not to force column vectors 
#'   to have a second dimension equal to one. Defaults to \code{FALSE}. This 
#'   behavior can also be changed globally using, for example 
#'   \code{options(atleast_2d = TRUE)}.
#' @return A  matrix or array of pseudorandom numbers. If the result is a 
#'         matrix, it will have class \code{c("matrix", "mat")}.
#' @seealso \code{\link{rand}}, \code{\link{randi}}, \code{\link{rnorm}}.
#' @export
#' @importFrom stats rnorm
#' @examples
#' randn(100, 100)  # 100 by 100 matrix of standard normal random variates
#' randn(2, 3, mean = 10, sd = 0.1)
randn <- function(nrow = 1, ncol = 1, ..., mean = 0, sd = 1, 
                  atleast_2d = NULL) {
  if (length(list(...)) == 0) {
    if (is.null(atleast_2d)) atleast_2d <- getOption("atleast_2d")
    if (!atleast_2d && ncol == 1) {  # no dim attribute!
      rnorm(nrow, mean = mean, sd = sd)
    } else {
      matrix(rnorm(nrow * ncol, mean = mean, sd = sd), nrow = nrow, 
             ncol = ncol)
    }
  } else {
    array(rnorm(nrow * ncol * prod(unlist(list(...))), mean = mean, sd = sd), 
          dim = c(nrow, ncol, unlist(list(...))))
  }
}


#' Resize Matrices and Arrays
#' 
#' Change shape and size of a matrix or array.
#' 
#' @param x A matrix or multi-way array.
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @param ... Further dimensions of the array.
#' @param across Character string specifying whether to flatten the matrix 
#'               across \code{"rows"} (default) or \code{"columns"}. This option
#'               is ignored for multi-way arrays.
#' @param byrow Logical. If FALSE (default) the new matrix is filled by columns, 
#'                       otherwise it is filled by rows. This option is ignored
#'                       for multi-way arrays.
#'              
#' @return A matrix of class \code{c("matrix", "mat")} with dimension 
#'         \code{nrow}-by-\code{ncol}.
#' @seealso \code{\link{flatten}}, \code{\link{mat}}, \code{\link{matrix}}.
#' @export
#' @examples
#' m <- 1:9
#' resize(m)
#' resize(m, 3, 3)
#' resize(m, 2, 2)
resize <- function(x, nrow, ncol, ..., across = c("rows", "columns"), 
                   byrow = FALSE) {
  
  ## Make sure x is a matrix. If it's not, try converting it.
  if (!is.matrix(x)) x <- as.matrix(x)
  if (missing(nrow)) nrow <- dim(x)[1L]  # keep first dimension
  if (missing(ncol)) ncol <- dim(x)[2L]  # keep second dimension
  
  # Flatten and reshape/resize matrix.
  across <- match.arg(across)
  if (length(list(...)) == 0) {
    x <- matrix(flatten(x, across = across), nrow = nrow, ncol = ncol, 
                byrow = byrow)
  } else {
    dim(x) <- c(nrow, ncol, unlist(list(...)))
  }
  x

}


#' Dimensions of a Matrix/Array
#' 
#' Retrieve the dimensions of a matrix or array.
#' 
#' @param x A matrix, array, or data frame.
#' @return The dimensions of the object.
#' @export
#' @seealso \code{\link{dim}}.
#' @examples
#' m <- mat("1, 3, 5; 7, 9, 11")
#' size(m)
size <- function(x) {
  dim(x)
}


#' Trace of a Matrix
#' 
#' Sum of diagonal elements of a matrix.
#' 
#' @param x A matrix.
#' @return The sum of the diagonal elements of \code{x}.
#' @export
#' @examples
#' tr(ones(5, 10))
#' x <- replicate(1000, tr(rand(25, 25)))
#' hist(x)
tr <- function(x) {
  sum(diag(x))  # sum of diagonal elements
}


#' Lower/Upper Triangular Matrix
#' 
#' Construct a matrix with ones at and below the given diagonal and zeros 
#' elsewhere.
#' 
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @param k The sub-diagonal at and below which the matrix is filled. 
#'          \code{k = 0} is the main diagonal, while k < 0 is below it, and 
#'          k > 0 is above. The default is 0.
#' @param diag Logical indicating whether to include the diagonal. Default is 
#'             \code{TRUE}.
#' @export
#' @examples
#' tri(5, 5)
#' tri(5, 5, 2)
#' tri(5, 5, -1)
tri <- function(nrow, ncol = nrow, k = 0, diag = TRUE) {
  # FIXME: Add an "upper = TRUE" option which would return t(m) instead?
  x <- matrix(nrow = nrow, ncol = ncol)
  if (diag) {
    m <- as.integer(row(x) >= col(x) - k)
  }
  else {
    m <- as.integer(row(x) > col(x) - k)
  }
  dim(m) <- dim(x)
  m
}


#' Extract Lower Triangular Matrix
#' 
#' Extract the lower triangular part of a matrix.
#' 
#' @param x A matrix.
#' @param k Diagonal above which to zero elements. \code{k = 0} (the default) is 
#'          the main diagonal, k < 0 is below it and k > 0 is above.
#' @param diag Logical indicating whether to include the diagonal. Default is 
#'             \code{TRUE}.
#' @export
#' @examples
#' tril(ones(5, 5))
#' tril(ones(5, 5), diag = TRUE)
tril <- function(x, k = 0, diag = TRUE) {
  if (diag) {
    m <- ifelse(row(x) >= col(x) - k, x, 0)
  }
  else {
    m <- ifelse(row(x) > col(x) - k, x, 0)
  }  
  class(m) <- typeof(x)
  dim(m) <- dim(x)
  m
}


#' Extract Upper Triangular Matrix
#' 
#' Extract the upper triangular part of a matrix.
#' 
#' @param x A matrix.
#' @param k Diagonal below which to zero elements. \code{k = 0} (the default) is 
#'          the main diagonal, k < 0 is below it and k > 0 is above.
#' @param diag Logical indicating whether to include the diagonal. Default is 
#'             \code{TRUE}.
#' @export
#' @examples
#' triu(ones(5, 5))
#' triu(ones(5, 5), diag = FALSE)
triu <- function(x, k = 0, diag = TRUE) {
  if (diag) {
    m <- ifelse(row(x) <= col(x) - k, x, 0)
  }
  else {
    m <- ifelse(row(x) < col(x) - k, x, 0)
  }  
  class(m) <- typeof(x)
  dim(m) <- dim(x)
  m
}


#' Lower Triangular Matrix Test
#' 
#' Determine if a Matrix is Lower Triangular
#' @param x A matrix
#' @return Logical indicating whether the given matrix is lower triangular.
#' @export
#' @examples
#' m <- mat("1, 0, 0, 0; -1, 1, 0, 0; -2, -2, 1, 0; -3, -3, -3, 1")
#' is.tril(m)
#' is.tril(eye(3, 5))
is.tril <- function(x) {
  # A matrix is lower triangular if all elements above the main diagonal are 
  # zero. Any number of the elements on the main diagonal can also be zero.
  all(x == tril(x))
}


#' Upper Triangular Matrix Test
#' 
#' Determine if a Matrix is Upper Triangular
#' @param x A matrix
#' @return Logical indicating whether the given matrix is lower triangular.
#' @export
#' @examples
#' m <- mat("1, -1, -1, -1; 0, 1, -2, -2; 0, 0, 1, -3; 0, 0, 0, 1")
#' is.triu(m)
#' is.triu(eye(3, 5))
is.triu <- function(x) {
  # A matrix is upper triangular if all elements below the main diagonal are 
  # zero. Any number of the elements on the main diagonal can also be zero.
  all(x == triu(x))
}