#' @title List with vectors
#' @description 
#' \code{list_with_vectors} checks if an object is a list with vectors \cr
#' \code{list_with_numeric_vectors} checks if an object is a 
#' list with numeric vectors \cr
#' \code{list_with_string_vectors} checks if an object is a 
#' list with string vectors
#' 
#' @param x an R object
#' @name list_with_vectors
#' @export list_with_vectors list_with_numeric_vectors
#' list_with_string_vectors
#' @aliases list_with_vectors list_with_numeric_vectors
#' list_with_string_vectors
#' @seealso \code{\link{is_vector}}, \code{\link{list_of_vectors}}
#' @examples
#' a = list(1:3, letters[1:3], c(exp(1), pi), NA)
#' b = list(1:3, c(exp(1), pi))
#' d = list(letters[1:3], 'bonjour a tous')
#' e = list(matrix(1:6, 2, 3), a, b)
#' 
#' list_with_vectors(1:10) # FALSE
#' list_with_vectors(b) # TRUE
#' list_with_vectors(d) # TRUE
#' 
#' list_with_numeric_vectors(a) # TRUE
#' list_with_numeric_vectors(b) # TRUE
#' list_with_string_vectors(d) # FALSE
#' 
#' list_with_string_vectors(a) # TRUE
#' list_with_string_vectors(d) # TRUE
#' list_with_string_vectors(b) # FALSE
NULL

list_with_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is.vector))
    if (sum(vectors) > 0) TRUE else FALSE    
  } else FALSE
}

list_with_numeric_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_numeric_vector))
    if (sum(vectors) > 0) TRUE else FALSE    
  } else FALSE
}

list_with_string_vectors <- function(x) {
  if (is.list(x)) {
    vectors = unlist(lapply(x, is_string_vector))
    if (sum(vectors) > 0) TRUE else FALSE    
  } else FALSE
}
