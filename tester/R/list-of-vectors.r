#' @title List of vectors
#' @description 
#' \code{list_of_vectors} checks if an object is a list of vectors \cr
#' \code{list_of_numeric_vectors} checks if an object is a 
#' list of numeric vectors \cr
#' \code{list_of_string_vectors} checks if an object is a 
#' list of string vectors
#' \code{list_of_logical_vectors} checks if an object is a 
#' list of logical vectors
#' 
#' @param x an R object
#' @name list_of_vectors
#' @export list_of_vectors list_of_numeric_vectors 
#' list_of_string_vectors list_of_logical_vectors
#' @aliases list_of_vectors list_of_numeric_vectors 
#' list_of_string_vectors list_of_logical_vectors
#' @seealso \code{\link{is_vector}}, \code{\link{list_with_vectors}}
#' @examples
#' a = list(1:3, letters[1:3], c(exp(1), pi), NA)
#' b = list(1:3, c(exp(1), pi))
#' d = list(letters[1:3], 'bonjour a tous')
#' e = list(matrix(1:6, 2, 3), a, b)
#' 
#' list_of_vectors(a) # TRUE
#' list_of_vectors(b) # TRUE
#' list_of_vectors(d) # TRUE
#' list_of_vectors(e) # FALSE
#' 
#' list_of_numeric_vectors(a) # FALSE
#' list_of_numeric_vectors(b) # TRUE
#' 
#' list_of_string_vectors(a) # FALSE
#' list_of_string_vectors(d) # TRUE
#' 
#' list_of_logical_vectors(a) # FALSE
#' list_of_logical_vectors(d) # TRUE
NULL

list_of_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is.vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}

list_of_numeric_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_numeric_vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}

list_of_string_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_string_vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}

list_of_logical_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_logical_vector))
    if (sum(vectors) == length(x)) TRUE else FALSE    
  } else FALSE
}
