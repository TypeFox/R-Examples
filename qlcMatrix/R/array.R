# =====
# Supplementary function for sparse tensors as Arrays
# extending the functionality of "spam" library
# =====

# use names in line with matrix/Matrix

sparseArray <- function(i, v = NULL, ... ) {
  if (is.null(v)) {
    simple_sparse_array(i, v = rep.int(1, times = nrow(i)), ...)
  } else {
    simple_sparse_array(i, v, ...)
  }
}

Array <- function(A) {
  if (is.data.frame(A)) {
    # assume 'long' format of factors, no values, so all values are "1"
    indices <- data.matrix(A)
    attr(indices, "dimnames") <- NULL
    sparseArray(i = indices
                , dimnames = sapply(A, levels)
                )
  } else {
    # assume A is an array
    as.simple_sparse_array(A)
  }
}

# turn "simple_triplet_matrix" from spam into "dgTMatrix"

as.Matrix <- function(M) {
  if (is.null(M$dimnames)) {
    n <- list(NULL, NULL)
  } else {
    n <- M$dimnames
  }
  sparseMatrix(i = M$i
               , j = M$j
               , x = M$v
               , dims = c(M$nrow, M$ncol)
               , dimnames = n
               , giveCsparse = FALSE
               )
}

# =====
# General function to unfold margins from sparse array
# unfolded margins get added as last margin of new array
#
# Speed is not optimized
# =====

unfold <- function(x, MARGINS) {

  ndim <- length(dim(x))

  if(max(MARGINS) > ndim) {
    stop("MARGINS larger than array size")
  }

  # keep unchanged dimensions
  old_coor <- x$i[ , -MARGINS, drop = FALSE]

  # make new coordinates and insert as final dimension
  # this is the sparse-magic, getting all coordinates right
  f <- head(cumprod(c(1, x$dim[MARGINS])), -1)
  make_new_coor <- function(coor) {
    1 + sum((coor - 1)*f)
  }

  new_coor <- apply(x$i[ , MARGINS, drop = FALSE], 1, make_new_coor)
  new_i <- cbind(old_coor, new_coor)

  # new size
  new_dim <- x$dim[ -MARGINS ]
  new_dim <- c(new_dim, prod(x$dim[MARGINS]))

  # permuation vector of position of new dimensions
  p <- 1:ndim
  p[-MARGINS] <- 1:(ndim-length(MARGINS))
  p[MARGINS] <- length(new_dim)

  # make new array
  a <- simple_sparse_array(i = new_i, v = x$v, dim = new_dim)
  attr(a, "permutation") <- p
  attr(a, "unfolded") <- MARGINS
  return(a)
}

# =====
# Special case of unfolding, result being a matrix
# this should emulate the "tenmat" function from Matlab Tensor Toolbox
# =====

unfold_to_matrix <- function(x, ROWS, COLS = NULL) {

  ndim <- length(dim(x))

  if (!is.null(COLS) && length(c(ROWS,COLS)) != ndim) {
    stop("ROWS and COLS must contain all margins of x")
  }

  if (is.null(COLS)) {
    COLS <- (1:ndim)[-ROWS]
  }
  
  if (length(ROWS) == 1) {
    unfoldC <- unfold(x, COLS)
  } else {
    unfoldR <- unfold(x, ROWS)
    unfoldC <- unfold(unfoldR, attr(unfoldR,"permutation")[COLS])
  }
  
  return(as.Matrix(as.simple_triplet_matrix(unfoldC)))
}

# ====
# for Matlab compatibility
# ====

tenmat <- unfold_to_matrix

