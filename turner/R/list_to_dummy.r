#' @title Create a dummy matrix from the elements in a list
#' 
#' @description
#' Create a dummy matrix based on the elements of a list. Each column in the
#' produced matrix is a dummy indicator.
#'
#' @param alist a list of vectors
#' @return A matrix of dummy variables
#' @author Gaston Sanchez
#' @seealso \code{\link{dummy_to_list}}, \code{\link{listify}}
#' @export
#' @examples
#' # let's say you have a list like this
#' num_list = list(1:3, 4:5, 6:9)
#' 
#' # get dummy matrix
#' list_to_dummy(num_list)
#' 
#' # try with a list of strings
#' str_list = list(c("a","b","c"), c("d", "e"), c("f","g","h"))
#' list_to_dummy(str_list)
list_to_dummy <- function(alist)
{
  if (!list_of_vectors(alist))
    stop("\n'list_to_dummy()' requires a list of vectors")
  
  aux = lengths(alist)
  to = cumsum(aux)
  from = to - aux + 1
  dummy_matrix = matrix(0, sum(aux), length(aux))
  for (j in seq_along(aux))
    dummy_matrix[from[j]:to[j], j] = 1
  dummy_matrix
}
