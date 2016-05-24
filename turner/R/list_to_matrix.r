#' @title Design-type matrix from the elements in a list
#' 
#' @description
#' Create a design-type matrix based on the elements of a list. Each column in
#' the produced matrix is linked to the vectors in the list. See example.
#'
#' @param alist a list of numeric vectors
#' @return A design-type matrix
#' @author Gaston Sanchez
#' @seealso \code{\link{list_to_dummy}}, \code{\link{indexify}}
#' @export
#' @examples
#' # let's say you have a list like this
#' num_list = list(1:3, 4:5, 6:9)
#' 
#' # get design-type matrix
#' list_to_matrix(num_list)
list_to_matrix <- function(alist)
{
  if (!list_of_numeric_vectors(alist))
    stop("\n'list_to_matrix()' requires a list of numeric vectors")
  
  aux = lengths(alist)
  to = cumsum(aux)
  from = to - aux + 1
  linked_matrix = matrix(0, sum(aux), length(aux))
  for (j in seq_along(aux))
    linked_matrix[from[j]:to[j], j] = alist[[j]]
  linked_matrix
}
