#' Set comparisons
#' 
#' Checks on the contents of two vectors (ignoring the order of the elements).
#' @param x A vector.
#' @param y Another vector.
#' @param strictly Logical.  If \code{TRUE}, \code{x} and \code{y} should not
#' be set equal.
#' @param .xname Not intended to be used directly.
#' @param .yname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{is_*} functions return \code{TRUE} or \code{FALSE}.
#' The \code{assert_*} functions throw an error in the event of failure.
#' @seealso \code{\link{is_subset}}, \code{\link[base]{sets}}, 
#' \code{\link[sets]{set_is_equal}}
#' @examples
#' # Same contents, different order, returns TRUE
#' are_set_equal(1:5, 5:1)
#' # Different lengths
#' are_set_equal(1:5, 1:6)
#' # First vector contains values not in second vector
#' are_set_equal(1:5, c(1:4, 4))
#' # Second vector contains values not in first vector
#' are_set_equal(c(1:4, 4), 1:5)
#' 
#' # Is x a subset of y?
#' is_subset(1:4, 1:5)
#' is_subset(1:5, 1:4)
#' 
#' # Is x a superset of y?
#' is_superset(1:5, 1:4)
#' is_superset(1:4, 1:5)
#' 
#' # The strictly argument checks for a strict sub/superset
#' is_subset(1:5, 1:5, strictly = TRUE)
#' is_superset(1:5, 1:5, strictly = TRUE)
#' 
#' # Types are coerced to be the same, as per base::setdiff
#' are_set_equal(1:4, c("4", "3", "2", "1"))
#' 
#' # Errors are thrown in the event of failure
#' assert_are_set_equal(1:5, 5:1)
#' assertive.base::dont_stop(assert_are_set_equal(1:5, 1:6))
#' 
#' assert_is_subset(1:4, 1:5)
#' assertive.base::dont_stop(assert_is_subset(1:5, 1:4))
#' 
#' assert_is_superset(1:5, 1:4)
#' assertive.base::dont_stop(assert_is_superset(1:4, 1:5))
#' 
#' # A common use case: checking that data contains required columns
#' required_cols <- c("Time", "weight", "Diet")
#' assert_is_superset(colnames(ChickWeight), required_cols)
#' @export
are_set_equal <- function(x, y, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y))
{
  x <- unique(x)
  y <- unique(y)
  if(length(x) != length(y))
  {
    return(
      false(
        gettext(
          "%s and %s have different numbers of elements (%d versus %d)."
        ), 
        .xname, 
        .yname,
        length(x),
        length(y)
      )
    )
  }
  if(!(ok <- is_subset(x, y, FALSE, .xname, .yname)))
  {
    return(ok)
  }  
  if(!(ok <- is_subset(y, x, FALSE, .yname, .xname)))
  {
    return(ok)
  }  
  TRUE
}

#' @rdname are_set_equal
#' @export
is_set_equal <- function(x, y, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y))
{
  .Deprecated("are_set_equal")
  are_set_equal(x, y, .xname, .yname)
}

#' @rdname are_set_equal
#' @export
is_subset <- function(x, y, strictly = FALSE, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y))
{
  diffxy <- setdiff(x, y)
  if(length(diffxy) > 0)
  {
    return(
      false(
        ngettext(
          length(diffxy), 
          "The element %s in %s is not in %s.", 
          "The elements %s in %s are not in %s."
        ),
        toString(sQuote(diffxy), 100),
        .xname,
        .yname
      )
    )
  } 
  if(strictly && length(setdiff(y, x)) == 0)
  {
    return(false("%s and %s are set equal.", .xname, .yname))
  }  
  TRUE
}

#' @rdname are_set_equal
#' @export
is_superset <- function(x, y, strictly = FALSE, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y))
{
  is_subset(y, x, strictly, .yname, .xname)
}
