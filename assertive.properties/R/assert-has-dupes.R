#' @rdname has_duplicates
#' @export
assert_has_duplicates <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                                
  assert_engine(
    has_duplicates, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname has_duplicates
#' @export
assert_has_no_duplicates <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                             
  assert_engine(
    has_no_duplicates, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
