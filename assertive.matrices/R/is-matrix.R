# TODO: checks for positive-definite (see 
# https://www.mathworks.com/matlabcentral/newsreader/view_thread/55350),
# unit, orthogonal matrices.

#' Is the input a diagonal matrix?
#' 
#' Checks that the input is a diagonal matrix.
#' 
#' @param x Input to check.
#' @param tol Absolute values smaller than \code{tol} are not considered.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if the input is all zeroes (after coercion to be a 
#' matrix).
#' @examples
#' x <- diag(3)
#' is_diagonal_matrix(x)
#' x[1, 2] <- 100 * .Machine$double.eps
#' is_diagonal_matrix(x)
#' x[2, 3] <- 101 * .Machine$double.eps
#' is_diagonal_matrix(x)
#' @export
is_diagonal_matrix <- function(x, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x))
{
  .xname <- force(.xname)
  x <- coerce_to(x, "matrix", .xname)
  diag(x) <- 0
  ok <- is_zero_matrix(x, tol = tol, .xname = .xname)
  if(!ok)
  {
    cause(ok) <- sub("non-zero", gettext("off-the-diagonal"), cause(ok))
    return(ok)
  }
  TRUE
}

#' Is the matrix an identity matrix?
#' 
#' Checks that the input is an identity matrix.
#' 
#' @param x Input to check.
#' @param tol Abolute deviations from the expected values smaller than 
#' \code{tol} are not considered.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if the input is all zeroes (after coercion to be a 
#' matrix).
#' @examples
#' x <- diag(3)
#' is_identity_matrix(x)
#' x[1, 2] <- 100 * .Machine$double.eps
#' is_identity_matrix(x)
#' x[2, 3] <- 101 * .Machine$double.eps
#' is_identity_matrix(x)
#' @export
is_identity_matrix <- function(x, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x))
{
  .xname <- force(.xname)
  x <- coerce_to(x, "matrix", .xname)
  if(!(ok <- is_square_matrix(x, .xname = .xname)))
  {
    return(ok)
  }
  diag(x) <- diag(x) - 1
  if(!(ok <- is_zero_matrix(x, tol = tol, .xname = .xname)))
  {
    cause(ok) <- sub(
      "%s contains", 
      gettext("%s, with 1 subtracted from the diagonal, contains"), 
      cause(ok)
    )
    return(ok)
  }
  TRUE
}

#' Is the matrix upper/lower triangular?
#' 
#' Checks that the input is an upper or lower triangular matrix.
#' 
#' @param x Input to check.
#' @param strictly Logical. If \code{TRUE}, the diagonal must consist of zeroes.
#' @param tol Abolute deviations from the expected values smaller than 
#' \code{tol} are not considered.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if the input is all zeroes (after coercion to be a 
#' matrix).
#' @examples
#' x <- matrix(c(1, 2, 3, 0, 5, 6, 0, 0, 9), nrow = 3)
#' is_lower_triangular_matrix(x)
#' is_lower_triangular_matrix(x, strictly = TRUE)
#' is_upper_triangular_matrix(t(x))
#' is_upper_triangular_matrix(t(x), strictly = TRUE)
#' x[1, 2] <- 100 * .Machine$double.eps
#' is_lower_triangular_matrix(x)
#' x[2, 3] <- 101 * .Machine$double.eps
#' is_lower_triangular_matrix(x)
#' @export
is_lower_triangular_matrix <- function(x, strictly = FALSE, 
  tol = 100 * .Machine$double.eps, .xname = get_name_in_parent(x))
{
  .xname <- force(.xname)
  x <- coerce_to(x, "matrix", .xname)
  x[lower.tri(x, diag = !strictly)] <- 0
  if(!(ok <- is_zero_matrix(x, tol = tol, .xname = .xname)))
  {
    cause(ok) <- paste(
      gettext("The upper triangular portion of"), 
      cause(ok)
    )
    return(ok)
  }
  TRUE
}

#' Is the matrix a square matrix?
#' 
#' Checks that the input is a square matrix.
#' 
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if the input is all zeroes (after coercion to be a 
#' matrix).
#' @examples
#' is_square_matrix(matrix(1:9, nrow = 3))
#' is_square_matrix(matrix(1:12, nrow = 3))
#' @export
is_square_matrix <- function(x, .xname = get_name_in_parent(x))
{
  .xname <- force(.xname)
  x <- coerce_to(x, "matrix", .xname)
  dimx <- dim(x)
  if(dimx[1L] != dimx[2L])
  {
    return(
      false(
        gettext("%s is not a square matrix; its dimensions are %s."),
        .xname, 
        toString(dimx)
      )
    )
  }
  TRUE
}

#' Is the input a symmetric matrix?
#'
#' Checks that the input is a symmetric matrix.
#' 
#' @param x Input to check.
#' @param tol Differences smaller than \code{tol} are not considered.
#' @param .xname Not intended to be used directly.
#' @param ... Passed to \code{all.equal}.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if the input is symmetric (after coercion to be a 
#' matrix).
#' @examples
#' m <- diag(3); m[3, 1] <- 1e-100
#' assert_is_symmetric_matrix(m)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_symmetric_matrix(m, tol = 0))
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base is_identical_to_true
#' @export
is_symmetric_matrix <- function(x, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x), ...)
{
  .xname <- force(.xname)
  x <- coerce_to(x, "matrix", .xname)
  if(!(ok <- is_square_matrix(x, .xname = .xname)))
  {
    return(ok)
  }
  symmetry_test <- if(is.complex(x)) 
  {
    all.equal.numeric(x, Conj(t(x)), tolerance = tol, ...)
  } else 
  {
    all.equal(x, t(x), tolerance = tol, ...)
  }
  if(!is_identical_to_true(symmetry_test))
  {
    return(false(gettext("%s is not a symmetric matrix."), .xname))
  }
  TRUE
}

#' @rdname is_lower_triangular_matrix
#' @export
is_upper_triangular_matrix <- function(x, strictly = FALSE, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x))
{
  .xname <- force(.xname)
  x <- coerce_to(x, "matrix", .xname)
  x[upper.tri(x, diag = !strictly)] <- 0
  if(!(ok <- is_zero_matrix(x, tol = tol, .xname = .xname)))
  {
    cause(ok) <- paste(
      gettext("The lower triangular portion of"), 
      cause(ok)
    )
    return(ok)
  }
  TRUE
}

#' Is the input a zero matrix
#' 
#' Checks that the input is a matrix of zeroes.
#' 
#' @param x Input to check.
#' @param tol Absolute values smaller than \code{tol} are not considered.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{TRUE} if the input is all zeroes (after coercion to be a 
#' matrix).
#' @examples
#' x <- matrix(numeric(9), 3)
#' is_zero_matrix(x)
#' x[1, 1] <- 100 * .Machine$double.eps
#' is_zero_matrix(x)
#' x[2, 2] <- 101 * .Machine$double.eps
#' is_zero_matrix(x)
#' @export
is_zero_matrix <- function(x, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x))
{
  .xname <- force(.xname)
  x <- coerce_to(x, "matrix", .xname)
  bad <- abs(x) > tol
  if(any(bad))
  {
    bad_data <- data.frame(which(bad, TRUE), value = x[bad])
    return(
      false(
        gettextf(
          "%s contains non-zero elements:\n%s", 
          .xname, 
          assertive.base:::print_and_capture(bad_data)
        )
      )
    )
  }
  TRUE
}
