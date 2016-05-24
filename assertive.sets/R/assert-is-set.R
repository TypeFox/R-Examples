#' @rdname are_set_equal
#' @export
assert_are_set_equal <- function(x, y, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    are_set_equal, 
    x, 
    y = y, 
    .xname = get_name_in_parent(x), 
    .yname = get_name_in_parent(y),
    severity = severity
  )
}

#' @rdname are_set_equal
#' @export
assert_is_set_equal <- function(x, y, 
  severity = getOption("assertive.severity", "stop"))
{
  .Deprecated("assert_are_set_equal")
  assert_are_set_equal(x, y, severity = severity)
}

#' @rdname are_set_equal
#' @export
assert_is_subset <- function(x, y, strictly = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    is_subset, 
    x, 
    y = y, 
    strictly = strictly, 
    .xname = get_name_in_parent(x), 
    .yname = get_name_in_parent(y),
    severity = severity
  ) 
}

#' @rdname are_set_equal
#' @export
assert_is_superset <- function(x, y, strictly = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    is_superset, 
    x, 
    y = y, 
    strictly = strictly, 
    .xname = get_name_in_parent(x), 
    .yname = get_name_in_parent(y),
    severity = severity
  ) 
}
