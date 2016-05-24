#' Are the inputs the same length/dimension?
#' 
#' Checks if the inputs are the same length, or have the same dimensions.
#' @param x An R object or expression.
#' @param y Another R object or expression.
#' @param ... Some R expressions.
#' @param l A list of R expressions.
#' @param .xname Not intended to be used directly.
#' @param .yname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{are_same_length} and \code{have_same_dims} return \code{TRUE}
#' if \code{x} and \code{y} have the same length, or their dimensions are 
#' identical.  The \code{assert_*} functions throw an error on failure.
#' 
#' The legacy function \code{are_same_length_legacy} allows an arbitrary number
#' of inputs and returns a symmetric square logical matrix which is \code{TRUE} 
#' where pairs of inputs are the same length.  (The new version of the function 
#' is easier to work with, and it is recommended that you switch your code to 
#' it.)
#' @seealso \code{\link[base]{length}}, \code{\link[assertive.base]{are_identical}}
#' @examples
#' are_same_length(runif(5), list(1, 2:3, 4:6, 7:10, 11:15))
#' assertive.base::dont_stop(
#'   assert_are_same_length(runif(6), list(1, 2:3, 4:6, 7:10, 11:15))
#' )
#' have_same_dims(
#'   matrix(1:12, nrow = 4), 
#'   data.frame(x = 1:4, y = 5:8, y = 9:12)
#' )
#' have_same_dims(1:5, matrix(1:5))
#' @export
are_same_length <- function(x, y, .xname = get_name_in_parent(x),
  .yname = get_name_in_parent(y))
{
  len_x <- length(x)
  len_y <- length(y)
  if(len_x != len_y)
  {
    return(
      false(
        gettext("%s has length %d but %s has length %d."),
        .xname,
        len_x,
        .yname,
        len_y
      )
    )
  }
  TRUE
}

get_dim_string <- function(x)
{
  if(is.null(x)) "NULL" else toString(x)
}

#' @rdname are_same_length
#' @export
have_same_dims <- function(x, y, .xname = get_name_in_parent(x),
  .yname = get_name_in_parent(y))
{
  dim_x <- dim(x)
  dim_y <- dim(y)
  if(!identical(dim_x, dim_y))
  {
    return(
      false(
        gettext("%s has dim %s but %s has dim %s."),
        .xname,
        get_dim_string(dim_x),
        .yname,
        get_dim_string(dim_y)
      )
    )
  }
  TRUE
}

#' @rdname are_same_length
#' @export
are_same_length_legacy <- function(..., l = list())
{
  .Deprecated("are_same_length")
  envir <- parent.frame()
  inputs <- as.list(match.call())[-1]
  inputs_in_list <- as.list(inputs$l)[-1]
  inputs <- c(inputs[names(inputs) != "l"], inputs_in_list)
  input_pairs <- expand.grid(expr1 = inputs, expr2 = inputs)
  equality <- apply(
    input_pairs, 
    1, 
    function(row)
    {       
      with(
        row,         
        length(eval(expr1, envir = envir)) == length(eval(expr2, envir = envir))
      )
    }
  )
  matrix(
    equality,
    nrow     = length(inputs),
    dimnames = list(inputs, inputs) 
  )
}
