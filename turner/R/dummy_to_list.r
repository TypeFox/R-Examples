#' @title Create an indexed list from a dummy matrix
#' 
#' @description
#' Create an indexed list from the columns of a dummy (or semi-dummy) matrix
#'
#' @param Dummy matrix (dummy by columns)
#' @return A list of indexed vectors
#' @author Gaston Sanchez
#' @seealso \code{\link{list_to_dummy}}, \code{\link{listify}}
#' @export
#' @examples
#' # let's say you have a list like this
#' some_list = list(1:3, 1:2, 1:4)
#' 
#' # first create a dummy matrix based on some_list
#' some_dummy = list_to_dummy(some_list)
#' 
#' # now apply 'dummy_to_list'
#' dummy_to_list(some_dummy)
#' 
#' # a semi-dummy matrix
#' semi_dummy = some_dummy
#' semi_dummy[semi_dummy != 0] = rnorm(listsize(some_list))
#' dummy_to_list(semi_dummy)
dummy_to_list <- function(Dummy)
{
  if (is_not_matrix(Dummy))
    stop("\n'dummy_to_list()' requires a (dummy) matrix")
  
  indices = apply(Dummy, 2, function(x) sum(x != 0))
  listify(indices)
}
