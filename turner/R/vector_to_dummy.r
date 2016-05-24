#' @title Create a dummy matrix from the elements in a vector
#' 
#' @description
#' Create a dummy matrix based on the elements of a vector. Each column in the
#' produced matrix is a dummy indicator.
#'
#' @param avector a numeric vector
#' @return A matrix of dummy variables
#' @author Gaston Sanchez
#' @seealso \code{\link{list_to_dummy}}, \code{\link{factor_to_dummy}}
#' @export
#' @examples
#' # let's say you have a list like this
#' num_vec = c(2, 3, 1, 4)
#' 
#' # get dummy matrix
#' vector_to_dummy(num_vec)
vector_to_dummy <- function(avector)
{
  if (!is_numeric_vector(avector))
    stop("\n'vector_to_dummy()' requires a numeric vector")
  
  num_rows = sum(avector)
  num_cols = length(avector)
  
  # starting-and-ending positions
  start_end = from_to(avector)
  from = start_end$from
  to = start_end$to
  
  # build dummy matrix
  dummy_matrix = matrix(0, num_rows, num_cols)
  for (k in 1:num_cols) {
    dummy_matrix[from[k]:to[k],k] = 1
  }
  dummy_matrix
}
