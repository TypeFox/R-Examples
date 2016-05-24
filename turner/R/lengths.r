#' @title Length of each element within a list
#' 
#' @description
#' Get the length of the elements contained in a list.
#'
#' @param alist a list
#' @param out string indicating the format of the output (\code{"vector"} or 
#' \code{"list"})
#' @return A vector (or list) with the lengths of the elements in \code{alist}
#' @author Gaston Sanchez
#' @seealso \code{\link{length}}, \code{\link{funlist}}
#' @export
#' @examples
#' # say you have some list
#' some_list = list(1:3, 4:5, 6:9)
#' 
#' # length of each vector (output in vector format)
#' lengths(some_list)
#' 
#' # length of each vector (output in list format)
#' lengths(some_list, out = 'list')
#' 
#' # compare to 'length()'
#' length(some_list)
lengths <- function(alist, out = "vector")
{
  if (!is.list(alist)) 
    stop("\n'lengths()' requires a list")
  
  bad_out <- !(out %in% c('vector', 'list'))
  if (bad_out) out = 'vector'
  
  if (out == "vector")
    unlist(lapply(alist, length))
  else lapply(alist, length)
}
