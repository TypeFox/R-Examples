#' How does the input relate to a value?
#'
#' Is \code{x} equal/not equal/greater than/less than \code{y}?
#' @param x A numeric vector.
#' @param y Another numeric vector, typically scalar or the same
#' length as \code{x}.  See note.
#' @param tol Values within \code{tol} are considered equal.
#' @param na_ignore A logical value.  If \code{FALSE}, \code{NA} values
#' cause an error; otherwise they do not.  Like \code{na.rm} in many
#' stats package functions, except that the position of the failing
#' values does not change.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @param .yname Not intended to be used directly.
#' @return \code{TRUE} if the input \code{x} is equal/not equal/greater 
#' than/less than \code{y}
#' @note The usual recycling rules apply when \code{x} and \code{y}
#' are different lengths.  See Intro to R for details:
#' \url{https://cran.r-project.org/doc/manuals/r-devel/R-intro.html#Vector-arithmetic}
#' \url{https://cran.r-project.org/doc/manuals/r-devel/R-intro.html#The-recycling-rule} 
#' @examples
#' # Approximate and exact floating point comparisons:
#' # See FAQ on R 7.31
#' x <- sqrt(2) * sqrt(2)
#' is_equal_to(x, 2)
#' is_equal_to(x, 2, tol = 0)
#' is_not_equal_to(x, 2)
#' is_not_equal_to(x, 2, tol = 0)
#' 
#' # Elements of x and y are recycled
#' is_equal_to(1:6, 1:3)
#' 
#' # Inequalities
#' x <- c(1 - .Machine$double.neg.eps, 1, 1 + .Machine$double.eps)
#' is_greater_than(x, 1)
#' is_greater_than_or_equal_to(x, 1)
#' is_less_than(x, 1)
#' is_less_than_or_equal_to(x, 1)
#' @export
is_equal_to <- function(x, y, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x), .yname = get_name_in_parent(x))
{
  if(!is.integer(x))
  {
    x <- coerce_to(x, "numeric", .xname)
  }
  if(!is.integer(y))
  {
    y <- coerce_to(y, "numeric", .yname)
  }
  difference <- abs(x - y)
  ok <- difference <= tol
  call_and_name(
    function(x) 
    {
      set_cause(ok, sprintf("not equal to %g (tol = %g); abs. difference = %g", y, tol, difference))
    }, 
    rep_len(x, length(ok))
  ) 
}

#' @rdname is_equal_to
#' @export
is_not_equal_to <- function(x, y, tol = 100 * .Machine$double.eps, 
  .xname = get_name_in_parent(x), .yname = get_name_in_parent(x))
{
  if(!is.integer(x))
  {
    x <- coerce_to(x, "numeric", .xname)
  }
  if(!is.integer(y))
  {
    y <- coerce_to(y, "numeric", .yname)
  }
  ok <- abs(x - y) > tol
  call_and_name(
    function(x) 
    {
      set_cause(ok, sprintf("equal to %g (tol = %g)", y, tol))
    }, 
    rep_len(x, length(ok))
  ) 
}

#' @rdname is_equal_to
#' @export
is_greater_than <- function(x, y, 
  .xname = get_name_in_parent(x), .yname = get_name_in_parent(x))
{
  if(!is.integer(x))
  {
    x <- coerce_to(x, "numeric", .xname)
  }
  if(!is.integer(y))
  {
    y <- coerce_to(y, "numeric", .yname)
  }
  ok <- x > y
  call_and_name(
    function(x) 
    {
      set_cause(ok, paste("less than or equal to", y))
    }, 
    rep_len(x, length(ok))
  ) 
}

#' @rdname is_equal_to
#' @export
is_greater_than_or_equal_to <- function(x, y, 
  .xname = get_name_in_parent(x), .yname = get_name_in_parent(x))
{
  if(!is.integer(x))
  {
    x <- coerce_to(x, "numeric", .xname)
  }
  if(!is.integer(y))
  {
    y <- coerce_to(y, "numeric", .yname)
  }
  ok <- x >= y
  call_and_name(
    function(x) 
    {
      set_cause(ok, paste("less than", y))
    }, 
    rep_len(x, length(ok))
  ) 
}

#' @rdname is_equal_to
#' @export
is_less_than <- function(x, y, 
  .xname = get_name_in_parent(x), .yname = get_name_in_parent(x))
{
  if(!is.integer(x))
  {
    x <- coerce_to(x, "numeric", .xname)
  }
  if(!is.integer(y))
  {
    y <- coerce_to(y, "numeric", .yname)
  }
  ok <- x < y
  call_and_name(
    function(x) 
    {
      set_cause(ok, paste("greater than or equal to", y))
    }, 
    rep_len(x, length(ok))
  ) 
}

#' @rdname is_equal_to
#' @export
is_less_than_or_equal_to <- function(x, y, 
  .xname = get_name_in_parent(x), .yname = get_name_in_parent(x))
{
  if(!is.integer(x))
  {
    x <- coerce_to(x, "numeric", .xname)
  }
  if(!is.integer(y))
  {
    y <- coerce_to(y, "numeric", .yname)
  }
  ok <- x <= y
  call_and_name(
    function(x) 
    {
      set_cause(ok, paste("greater than", y))
    }, 
    rep_len(x, length(ok))
  ) 
}
