#' @title Size: total number of elements in a list
#' 
#' @description
#' Get the total number of elements in a list.
#'
#' @param alist a list
#' @return number of elements in \code{alist}.
#' @author Gaston Sanchez
#' @seealso \code{\link{lengths}}
#' @aliases listsize sizelist
#' @export listsize sizelist
#' @examples
#' some_list = list(1:3, 4:5, 6:9)
#' 
#' # number of elems in 'some_list'
#' listsize(some_list)
listsize <- sizelist <- function(alist)
{
  if (!is.list(alist))
    stop("\n'listsize()' requires a list")
  
  length(unlist(alist))
}
