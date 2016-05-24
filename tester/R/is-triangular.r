#' @title Is triangular matrix
#' 
#' @description
#' \code{is_lower_triangular} tests if a matrix is lower triangular \cr
#' \code{is_upper_triangular} tests if a matrix is upper triangular \cr
#' \code{is_triangular_matrix} tests if a matrix is triangular (both
#' lower or upper triangular)
#' 
#' @param x a matrix
#' @param diag should the diagonal be included? (\code{FALSE} by default)
#' @name is_triangular_matrix
#' @aliases is_lower_triangular is_upper_triangular is_triangular_matrix
#' @export is_lower_triangular is_upper_triangular is_triangular_matrix
#' @examples
#' some_matrix = matrix(1:9, 3, 3)
#' lower_matrix <- upper_matrix <- some_matrix
#' lower_matrix[upper.tri(some_matrix)] <- 0
#' upper_matrix[lower.tri(some_matrix)] <- 0
#' 
#' is_triangular_matrix(some_matrix) # TRUE
#' is_triangular_matrix(lower_matrix) # TRUE
#' is_triangular_matrix(upper_matrix) # TRUE
#' 
#' is_lower_triangular(some_matrix) # FALSE
#' is_lower_triangular(lower_matrix) # FALSE
#' is_lower_triangular(upper_matrix) # FALSE
#' 
#' is_upper_triangular(some_matrix) # FALSE
#' is_upper_triangular(lower_matrix) # FALSE
#' is_upper_triangular(upper_matrix) # FALSE
NULL

is_lower_triangular <- function(x, diag = FALSE) {
  if (is.matrix(x)) {
    all(x[upper.tri(x, diag = diag)] == 0)
  } else FALSE
}

is_upper_triangular <- function(x, diag = FALSE) {
  if (is.matrix(x)) {
    all(x[lower.tri(x, diag = diag)] == 0)
  } else FALSE
}

is_triangular_matrix <- function(x, diag = FALSE) {
  if (is.matrix(x)) {
    is_lower_triangular(x) | is_upper_triangular(x)
  } else FALSE
}
