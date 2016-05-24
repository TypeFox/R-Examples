#' Does the input have rows/columns?
#'
#' Checks to see if the input has rows/columns.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_rows} and \code{has_cols} return \code{TRUE} if 
#' \code{nrow} and \code{ncol} respectively return a value that is 
#' non-null and positive.  The \code{assert_*} functions return nothing 
#' but throw an error if the corresponding \code{has_*} function returns
#' \code{FALSE}.
#' @seealso \code{\link{ncol}}.
#' @examples
#' assert_has_rows(data.frame(x = 1:10))
#' assert_has_cols(matrix())
#' @export
has_cols <- function(x, .xname = get_name_in_parent(x))
{
  ncolx <- ncol(x)
  if(is.null(ncolx)) 
  {
    return(false("The number of columns in %s is NULL.", .xname))  
  }
  if(ncolx == 0L) 
  {
    return(false("The number of columns in %s is zero.", .xname))
  }
  TRUE
} 

#' Does the input have dimensions?
#'
#' Checks to see if the input has dimensions.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_dims} returns\code{TRUE} if \code{dim} is non-null.
#' \code{assert_has_dims} returns nothing but throws an error if
#' \code{has_dims} is not \code{TRUE}.
#' @seealso \code{\link[base]{dim}}, \code{\link{is_of_dimension}}.
#' @export
has_dims <- function(x, .xname = get_name_in_parent(x))
{
  dim_x <- dim(x)
  if(is.null(dim_x)) 
  {
    return(false("The dimensions of %s are NULL.", .xname))
  }
  TRUE
}

#' @rdname has_cols
#' @export
has_rows <- function(x, .xname = get_name_in_parent(x))
{
  nrowx <- nrow(x)
  if(is.null(nrowx)) 
  {
    return(false("The number of rows in %s is NULL.", .xname))  
  }
  if(nrowx == 0L) 
  {
    return(false("The number of rows in %s is zero.", .xname))
  }
  TRUE
} 
