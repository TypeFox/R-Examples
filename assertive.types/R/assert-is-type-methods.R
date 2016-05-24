#' @rdname is_class
#' @export
assert_all_are_classes <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_class, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_class
#' @export
assert_any_are_classes <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_class, 
    x, 
    .xname = get_name_in_parent(x),
    what = "any",
    severity = severity
  )
}
