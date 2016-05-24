###########################################################################/**
# @RdocFunction snpArrayToMatrix
# @alias snpMatrixToArray
#
# @title "Reshapes SNP data in matrix form to array form and vice versa"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{Y}{A 2LxI @matrix or a Lx2xI @array,
#    where L is the number of probe pairs and I is the number of arrays.}
#  \item{dropNames}{If @TRUE, dimension names are dropped,
#    otherwise preserved.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a Lx2xI @array or a 2LxI matrix.
# }
#
# @examples "../incl/snpArrayToMatrix.Rex"
#
# @keyword internal
#*/###########################################################################
snpArrayToMatrix <- function(Y, dropNames=TRUE, ...) {
  # Argument 'Y':
  if (!is.array(Y)) {
    stop("Argument 'Y' is not an array: ", class(Y)[1L]);
  }
  dim <- dim(Y);
  if (length(dim) != 3L) {
    dimStr <- paste(dim, collapse="x");
    stop("Argument 'Y' is not a three-dimensional array: ", dimStr);
  }

  if (dim[2] != 2L) {
    stop("Second dimension of argument 'Y' is not of length 2: ", dim[2L]);
  }

  # Transform to a "stacked" matrix
  dim <- c(dim[1L]*dim[2L], dim[3L]);
  if (dropNames) {
    dim(Y) <- dim;
  } else {
    dimnames <- list(NULL, dimnames(Y)[[3L]]);
    dim(Y) <- dim;
    dimnames(Y) <- dimnames;
  }

  Y;
} # snpArrayToMatrix()


snpMatrixToArray <- function(Y, dropNames=TRUE, ...) {
  # Argument 'Y':
  if (!is.matrix(Y)) {
    stop("Argument 'Y' is not a matrix: ", class(Y)[1]);
  }
  dim <- dim(Y);
  if (dim[1L] %% 2 != 0L) {
    stop("The length of the first dimension of argument 'Y' is not even: ", dim[1]);
  }

  # Unstack to an array
  dim <- c(dim[1L]/2, 2L, dim[2L]);
  if (dropNames) {
    dim(Y) <- dim;
  } else {
    dimnames <- list(NULL, c("A", "B"), dimnames(Y)[[2L]]);
    dim(Y) <- dim;
    dimnames(Y) <- dimnames;
  }

  Y;
} # snpMatrixToArray()


############################################################################
# HISTORY:
# 2009-03-25 [HB]
# o Created.
############################################################################
