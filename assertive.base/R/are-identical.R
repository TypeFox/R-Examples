#' Are the inputs identical?
#' 
#' Checks if the inputs are identical.
#' @param x An R object or expression.
#' @param y Another R object or expression.
#' @param allow_attributes If \code{TRUE}, The attributes of \code{x} and 
#' \code{y} are allowed to differ.
#' @param ... Some R expressions, deprecated.
#' @param l A list of R expressions, deprecated.
#' @param .xname Not intended to be used directly.
#' @param .yname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{are_identical} returns \code{TRUE} if \code{x} and \code{y} 
#' are identical.  The \code{assert_*} function throws an error on failure.
#' 
#' The legacy function \code{are_identical_legacy} allows an arbitrary number
#' of inputs and returns a symmetric square logical matrix which is \code{TRUE} 
#' where pairs of inputs are identical.  (The new version of the function 
#' is easier to work with, and it is recommended that you switch your code to 
#' it.)
#' @seealso \code{\link[base]{identical}}, 
#' \code{\link[assertive.properties]{are_same_length}}
#' @examples
#' x <- 1:5
#' are_identical(c(1, -1), cos(c(0, pi)))
#' assertive.base::dont_stop(assert_are_identical(c(1, 1), cos(c(0, pi))))
#' @export
are_identical <- function(x, y, allow_attributes =FALSE,
  .xname = get_name_in_parent(x), .yname = get_name_in_parent(y))
{  
  if(allow_attributes) 
  {
    x <- strip_attributes(x)
    y <- strip_attributes(y)
  }
  if(!identical(x, y))
  {
    return(
      false(
        gettext("%s and %s are not identical."),
        .xname,
        .yname
      )
    )
  }
  TRUE
}



#' @rdname are_identical
#' @export
are_identical_legacy <- function(..., l = list())
{
  envir <- parent.frame()
  inputs <- as.list(match.call())[-1]
  inputs_in_list <- as.list(inputs$l)[-1]
  inputs <- c(inputs[names(inputs) != "l"], inputs_in_list)
  input_pairs <- expand.grid(expr1 = inputs, expr2 = inputs)
  identicality <- apply(
    input_pairs, 
    1, 
    function(row)
    {       
      with(
        row, 
        identical(
          eval(expr1, envir = envir),
          eval(expr2, envir = envir)
        )
      )
    }
  )
  matrix(
    identicality,
    nrow     = nargs(),
    dimnames = list(inputs, inputs) 
  )
}
