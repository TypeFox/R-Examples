#' @title Create indices for elements in a vector or list
#' 
#' @description
#' Create indexed components for the elements of a list.
#'
#' @param x a numeric vector or list of vectors
#' @param out string indicating the output format 
#' (\code{"vector"} or \code{"list"})
#' @return A vector (or list) of indexed numbers
#' @author Gaston Sanchez
#' @seealso \code{\link{listify}}
#' @export
#' @examples
#' # let's say you have a numeric vector like this
#' num_vec = c(2, 3, 1, 4)
#' 
#' # get indices in vector format
#' indexify(num_vec)
#' 
#' # let's say you have a list like this
#' str_list = list(c("a","b","c"), c("d", "e"), c("f","g","h"))
#' 
#' # get indices in vector format
#' indexify(str_list)
#' 
#' # get indices in list format
#' indexify(str_list, "list")
indexify <- function(x, out) {
  UseMethod("indexify", x)  
}


#' @S3method indexify default
indexify.default <- function(x, ...)
{
  if (!is_numeric_vector(x) || !list_of_vectors(x))
    stop("\n'indexify()' requires a numeric vector or a list of vectors")
}


#' @S3method indexify numeric
indexify.numeric <- function(x, out = "vector")
{
  if (!is_numeric_vector(x))
    stop("\n'indexify()' requires a numeric vector")
  
  if (out == "vector")
    rep(seq_along(x), x)
  else mapply(rep, seq_along(x), x)
}


#' @S3method indexify list
indexify.list <- function(x, out = "vector")
{
  if (!list_of_vectors(x))
    stop("\n'indexify()' requires a list of vectors")
  
  aux = unlist(lapply(x, length))
  if (out == "vector")
    rep(seq_along(aux), aux)
  else mapply(rep, seq_along(aux), aux)
}
