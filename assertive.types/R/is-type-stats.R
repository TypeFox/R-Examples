#' Is the input a (dendrogram) leaf?
#'
#' Checks to see if the input is a (dendrogram) leaf.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_leaf} reimplements \code{is.leaf}, providing more 
#' information on failure.
#' @seealso \code{\link[stats]{dendrogram}}.
#' @importFrom assertive.base is_identical_to_true
#' @export
is_leaf <- function(x, .xname = get_name_in_parent(x))
{
  leaf <- attr(x, "leaf")
  if(is.null(leaf)) 
  {
    return(false(gettext("%s has no 'leaf' attribute."), .xname))
  }
  ok <- is_identical_to_true(
    leaf, 
    allow_attributes = TRUE,
    paste("The leaf attribute of", .xname)
  )
  if(!ok)
  {
    return(ok)
  }
  TRUE
}

#' @rdname is_ts
#' @export
is_mts <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "mts", .xname)
}

#' @rdname is_function
#' @importFrom stats is.stepfun
#' @export
is_stepfun <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_function(x))) return(ok)
  if(!is.stepfun(x))
  {
    return(false(gettext("%s is not a step function."), .xname))
  }
  TRUE
} 

#' Is the input a time series?
#'
#' Checks to see if the input is a time series.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_ts} wraps \code{is.ts}, providing more 
#' information on failure.  \code{assert_is_ts} returns nothing but
#' throws an error if \code{is_ts} returns \code{FALSE}.
#' @seealso \code{\link[stats]{is.ts}}.
#' @examples
#' assert_is_ts(ts(1:10))
#' @importFrom assertive.base is2
#' @importFrom assertive.properties is_non_empty
#' @export
is_ts <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is2(x, "ts", .xname))) return(ok)
  if(!(ok <- is_non_empty(x, .xname = .xname))) return(ok)
  TRUE
}

#' @rdname is_ts
#' @importFrom assertive.base is2
#' @export
is_tskernel <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "tskernel", .xname)
}

