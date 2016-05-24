#' Is the input empty/scalar?
#'
#' Checks to see if the input has length zero/one.
#'
#' @param x Input to check.
#' @param n Non-negative integer(s) of the expected length/number of elements/
#' lengths of dimensions.  See note.
#' @param metric A string. Should be length or the number of elements be used to
#' determine if the object is empty/non-empty/scalar?
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_empty} returns \code{TRUE} if the input has length 
#' zero.  \code{is_scalar} returns \code{TRUE} if the input has length 
#' one.  The \code{assert_*} functions return nothing but throw an
#' error if the corresponding \code{is_*} function returns \code{FALSE}.
#' @note For \code{is_empty}, \code{is_non_empty} and \code{is_scalar}, \code{n}
#' should be an single integer representing either the expected length or the
#' expected number of elements in \code{x}.  For \code{is_of_dimension} \code{n}
#' should be a vector of integers representing the expected lengths of 
#' dimensions.
#' @seealso \code{\link{length}}.
#' @examples
#' # is_of_length returns TRUE if the length of an object
#' # matches a specified value.
#' is_of_length(1:5, 5)
#' assert_is_of_length(1:5, 5)
#' 
#' # has_elements returns TRUE if an object has a specified
#' # number of elements.  This is usually the same thing.
#' has_elements(1:5, 5)
#' assert_has_elements(1:5, 5)
#' 
#' # Data frames and lists behave differently for length
#' # and number of elements.
#' d <- data.frame(x = 1:5, y = letters[1:5])
#' assert_is_of_length(d, 2)
#' assert_has_elements(d, 10)
#' 
#' l <- list(a = 1:5, b = list(b.a = 1:3, b.b = 1:7))
#' assert_is_of_length(l, 2)
#' assert_has_elements(l, 15)
#' 
#' # Functions always have length one, but may have lots of 
#' # elements.
#' assert_is_of_length(var, 1)
#' assert_has_elements(var, 54)
#' 
#' # is_scalar is a shortcut for length one, or one elements.
#' assert_is_scalar(99)
#' assert_is_scalar("Multiple words in a single string are scalar.")
#' assert_is_scalar(NA)
#' 
#' # The two metrics can yield different results!
#' is_scalar(list(1:5))
#' is_scalar(list(1:5), "elements")
#' is_scalar(var)
#' is_scalar(var, "elements")
#' 
#' # Similarly, is_empty is a shortcut for length zero/zero elements.
#' assert_is_empty(NULL)
#' assert_is_empty(numeric())
#' assert_is_non_empty(1:10)
#' assert_is_non_empty(NA)
#' 
#' # is_of_dimension tests the lengths of all dimensions.
#' assert_is_of_dimension(d, c(5, 2))
#' assert_is_of_dimension(l, NULL)
#' @export
is_empty <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{  
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  metric_fn(x, 0L, .xname)
}


#' @rdname is_empty
#' @export
is_non_empty <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  if(metric_fn(x, 0)) 
  {
    msg <- switch(
      metric,
      length = gettext("%s has length 0."),
      elements = gettext("%s has 0 elements.")
    )
    return(false(msg, .xname))
  }
  TRUE
}

#' @rdname is_empty
#' @export
is_non_scalar <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  if(metric_fn(x, 1)) 
  {
    msg <- switch(
      metric,
      length = gettext("%s has length 1."),
      elements = gettext("%s has 1 element.")
    )
    return(false(msg, .xname))
  }
  TRUE
}

#' @rdname is_empty
#' @export
is_scalar <- function(x, metric = c("length", "elements"), 
                      .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  metric_fn(x, 1L, .xname)
}     


#' @rdname is_empty
#' @export
has_elements <- function(x, n, .xname = get_name_in_parent(x))
{
  n <- use_first(n)
  check_n(n)
  n_elements_x <- n_elements(x)
  if(n_elements_x != n)
  {
    return(
      false(
        ngettext(
          n_elements_x, 
          "%s has %d element, not %d.", 
          "%s has %d elements, not %d."
        ),
        .xname, 
        n_elements_x,
        n
      )
    )
  }
  TRUE
}

#' @rdname is_empty
#' @export
is_of_dimension <- function(x, n, .xname = get_name_in_parent(x))
{
  dim_x <- dim(x)
  # There are two cases to test: n is NULL, or n is a vector of natural 
  # numbers.
  if(is.null(n))
  {
    if(has_dims(x))
    {
      return(
        false(
          ngettext(
            length(dim_x), 
            "%s has dimension %s, not NULL.", 
            "%s has dimensions %s, not NULL."
          ), 
          .xname,
          deparse(dim_x)
        )
      )
    }
    return(TRUE)
  }
  check_n(n)
  if(!is_of_length(dim_x, length(n)))
  {
    return(
      false(
        ngettext(
          length(dim_x), 
          "%s has %d dimension, not %d.", 
          "%s has %d dimensions, not %d."
        ),  
        .xname, 
        length(dim_x),
        length(n)
      )
    )
  }
  differences <- dim_x != n
  if(any(differences))
  {
    bad <- which(differences)
    return(
      false(
        ngettext(
          length(bad), 
          "Dimension %s of %s is incorrect.", 
          "Dimensions %s of %s are incorrect."
        ), 
        toString(bad), 
        .xname
      )
    )
  }
  TRUE
}

#' @rdname is_empty
#' @export
is_of_length <- function(x, n, .xname = get_name_in_parent(x))
{
  n <- use_first(n)
  check_n(n)
  length_x <- length(x)
  if(length_x != n)
  {
    return(false("%s has length %d, not %d.", .xname, length_x, n))
  }
  TRUE
}

check_n <- function(n)
{
  if(n < 0 || n != round(n))
  {
    stop("n should be a non-negative integer vector.")
  }
}

get_metric <- function(metric)
{
  switch(
    metric,
    length   = is_of_length,
    elements = has_elements,
    stop("Bug in assertive; the metric", metric, "is not valid.", domain = NA)
  )
}
