#' @title Create a dummy matrix from the elements in a factor
#' 
#' @description
#' Create a dummy matrix based on the elements of a factor. Each column in the
#' produced matrix is a dummy indicator.
#'
#' @param afactor a factor (preferably of vectors)
#' @return A matrix of dummy variables
#' @author Gaston Sanchez
#' @seealso \code{\link{vector_to_dummy}}, \code{\link{list_to_dummy}}
#' @export
#' @examples
#' # let's say you have a list like this
#' some_factor = iris$Species[c(1:3,51:53,101:103)]
#' 
#' # get dummy matrix
#' factor_to_dummy(some_factor)
factor_to_dummy <- function(afactor)
{
  if (!is.factor(afactor))
    stop("\n'factor_to_dummy()' requires a factor")
 
  num_obs = length(afactor)
  categs = levels(afactor)
  num_categs = length(categs)
  obs_per_categ = tabulate(afactor)
  
  # build dummy matrix
  dummy_matrix = matrix(0, num_obs, num_categs)
  for (k in 1:num_categs) {
    tmp <- afactor == categs[k]
    dummy_matrix[tmp,k] = 1
  }
  colnames(dummy_matrix) = levels(afactor)
  rownames(dummy_matrix) = 1:num_obs
  # output
  dummy_matrix
}
