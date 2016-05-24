#' Get the dimensions of an object
#' 
#' Get the dimensions of an object, retuning the length if that object has no
#' \code{dim} attribute.
#' @param x Any object.
#' @return A integer vector of non-negative values.
#' @seealso \code{\link[base]{NROW}}, \code{\link[base]{dim}}
#' @examples
#' # For data frames and matrices, DIM is the same as dim.
#' DIM(sleep) 
#' # For vectors (and other objects without a dim attribute), DIM is the 
#' # same as length.
#' DIM(1:10)
#' DIM(list(x = 1:10))
#' @export
DIM <- function(x)
{
  dim_x <- dim(x)
  if(is.null(dim_x)) length(x) else dim_x
}

#' Get the number of elements
#' 
#' Gets the number of elements in an object.
#' @param x Any object.
#' @return A non-negative integer of the number of elements.
#' @note For atomic objects, the number of elements is the product of the
#' dimensions, as calculated by \code{\link{DIM}}.  For recursive objects,
#' the number of elements is the sum of the number of elements of each of
#' their atomic components.
#' @seealso \code{\link{DIM}}
#' @examples
#' n_elements(1:10)
#' n_elements(NULL)
#' n_elements(data.frame(x = 1:5, y = rnorm(5)))
#' n_elements(list(1:5, list(1:3, list(1:7))))
#' n_elements(var) # depends upon the length of the body
#' @export
n_elements <- function(x)
{
  if(is.recursive(x))
  {
    sum(vapply(x, n_elements, integer(1)))
  } else
  {
    as.integer(prod(DIM(x)))
  }  
}
