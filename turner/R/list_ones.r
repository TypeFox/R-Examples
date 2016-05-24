#' @title List with vectors of ones
#' 
#' @description
#' Create a list with vectors of ones from a numeric vector
#'
#' @param x a numeric vector
#' @return A list of vectors with ones
#' @author Gaston Sanchez
#' @seealso \code{\link{listify}}
#' @export
#' @examples
#' # let's say you have a numeric vector like this
#' num_vec = c(2, 3, 1, 4)
#' 
#' # get indices in vector format
#' list_ones(num_vec)
list_ones <- function(x)
{
  if (!is_numeric_vector(x))
    stop("\n'list_ones()' requires a numeric vector")
  # output
  lapply(x, function(u) rep(1, u))
}
