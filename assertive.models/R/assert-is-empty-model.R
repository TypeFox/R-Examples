#' @rdname is_empty_model
#' @export
assert_is_empty_model <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  assert_engine(
    is_empty_model, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_empty_model
#' @export
assert_is_non_empty_model <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  assert_engine(
    is_non_empty_model, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
