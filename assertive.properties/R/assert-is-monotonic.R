#' @rdname is_monotonic_increasing
#' @export
assert_is_monotonic_increasing <- function(x, strictly = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                    
  assert_engine(
    is_monotonic_increasing, 
    x, 
    strictly = strictly,
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_monotonic_increasing
#' @export
assert_is_monotonic_decreasing <- function(x, strictly = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                    
  assert_engine(
    is_monotonic_decreasing, 
    x, 
    strictly = strictly,
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
