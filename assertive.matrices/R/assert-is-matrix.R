#' @include imports.R

#' @rdname is_diagonal_matrix
#' @export
assert_is_diagonal_matrix <- function(x, tol = 100 * .Machine$double.eps,
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_diagonal_matrix, 
    x, 
    tol = tol, 
    .xname = get_name_in_parent(x),
    severity = severity
  )       
}

#' @rdname is_identity_matrix
#' @export
assert_is_identity_matrix <- function(x, tol = 100 * .Machine$double.eps, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_identity_matrix, 
    x, 
    tol = tol, 
    .xname = get_name_in_parent(x),
    severity = severity
  )       
}

#' @rdname is_lower_triangular_matrix
#' @export
assert_is_lower_triangular_matrix <- function(x, strictly = FALSE, tol = 100 * .Machine$double.eps, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_lower_triangular_matrix, 
    x,
    strictly = strictly,
    tol = tol, 
    .xname = get_name_in_parent(x),
    severity = severity
  )       
}

#' @rdname is_square_matrix
#' @export
assert_is_square_matrix <- function(x,
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_square_matrix, 
    x,
    .xname = get_name_in_parent(x),
    severity = severity
  )       
}

#' @rdname is_symmetric_matrix
#' @export
assert_is_symmetric_matrix <- function(x, tol = 100 * .Machine$double.eps, ..., 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_symmetric_matrix, 
    x, 
    tol = tol, 
    .xname = get_name_in_parent(x),
    ...
  )       
}

#' @rdname is_lower_triangular_matrix
#' @export
assert_is_upper_triangular_matrix <- function(x, strictly = FALSE, tol = 100 * .Machine$double.eps, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_upper_triangular_matrix, 
    x,
    strictly = strictly,
    tol = tol, 
    .xname = get_name_in_parent(x),
    severity = severity
  )       
}

#' @rdname is_zero_matrix
#' @export
assert_is_zero_matrix <- function(x, tol = 100 * .Machine$double.eps, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_zero_matrix, 
    x, 
    tol = tol, 
    .xname = get_name_in_parent(x),
    severity = severity
  )       
}
