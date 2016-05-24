#' @export
as.list.deque <- function(x, ...)
{
  .Call("R_deque_to_Rlist", x)
}



#' @export
head.deque <- function(x, n=6, ...)
{
  n <- as.integer(n)
  
  .Call("R_deque_headsortails", x, n, 1L)
  invisible()
}



#' @export
tail.deque <- function(x, n=6, ...)
{
  n <- as.integer(n)
  
  .Call("R_deque_headsortails", x, n, 2L)
  invisible()
}




#' Convert to Deque
#' 
#' @param x
#' An object either to be converted to the first element of a deque
#' (default), or the elements of a list (or columns of a dataframe)
#' to be set as elements of a deque.
#'
#' @return
#' A deque object.
#' 
#' @examples
#' library(dequer)
#' d <- as.deque(lapply(1:5, identity))
#' d
#' 
#' @export
#' @name as.deque
#' @rdname as.deque
as.deque <- function(x) UseMethod("as.deque")

#' @export
#' @rdname as.deque
as.deque.list <- function(x)
{
  d <- deque()
  
  for (obj in x)
    pushback(d, obj)
  
  return(d)
}

#' @export
#' @rdname as.deque
as.deque.default <- function(x)
{
  d <- deque()
  pushback(d, x)
  
  return(d)
}

