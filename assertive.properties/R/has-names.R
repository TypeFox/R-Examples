#' @rdname has_names
#' @export
has_colnames <- function(x, .xname = get_name_in_parent(x))
{
  colnamesx <- colnames(x)
  if(is.null(colnamesx)) 
  {
    return(false("The column names of %s are NULL.", .xname))
  }
  if(!any(nzchar(colnamesx))) 
  {
    return(false("The column names of %s are all empty.", .xname))
  }
  TRUE
} 

#' @rdname has_names
#' @export
has_dimnames <- function(x, .xname = get_name_in_parent(x))
{
  dimnamesx <- dimnames(x)
  if(is.null(dimnamesx)) 
  {
    return(false("The dimension names of %s are NULL.", .xname))
  }
  if(!any(nzchar(unlist(dimnamesx, use.names = FALSE)))) 
  {
    return(false("The dimension names of %s are all empty.", .xname))
  }
  TRUE
} 

#' Does the input have names?
#'
#' Checks to see if the input has (row/column/dimension) names.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_names} returns \code{TRUE} if \code{names} is 
#' non-null. 
#' \code{has_rownames}, \code{has_colnames} and \code{has_dimnames} work
#' in a similar fashion, checking the corresponding attributes.
#' \code{assert_has_names} returns nothing but throws an error if 
#' \code{has_names} is not \code{TRUE}.
#' @note Empty names (i.e., \code{""}) are not allowed in R, and are 
#' not checked here.
#' @seealso \code{\link[base]{names}}, \code{\link[base]{rownames}}, 
#' \code{\link[base]{colnames}}, \code{\link[base]{dimnames}}.
#' @examples
#' assert_has_names(c(a = 1, 2))
#' dfr <- data.frame(x = 1:5)
#' assert_has_rownames(dfr)
#' assert_has_colnames(dfr)
#' assert_has_dimnames(dfr)
#' @export
has_names <- function(x, .xname = get_name_in_parent(x))
{
  namesx <- names(x)
  if(is.null(namesx)) 
  {
    return(false("The names of %s are NULL.", .xname))
  }
  if(!any(nzchar(namesx))) 
  {
    return(false("The names of %s are all empty.", .xname))
  }
  TRUE
} 

#' @rdname has_names
#' @export
has_rownames <- function(x, .xname = get_name_in_parent(x))
{
  rownamesx <- rownames(x)
  if(is.null(rownamesx)) 
  {
    return(false("The row names of %s are NULL.", .xname))
  }
  if(!any(nzchar(rownamesx))) 
  {
    return(false("The row names of %s are all empty.", .xname))
  }
  TRUE
} 

