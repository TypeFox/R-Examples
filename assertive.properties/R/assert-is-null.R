#' @rdname is_null
#' @export
assert_is_not_null <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                      
  assert_engine(is_not_null, x, .xname = get_name_in_parent(x))   
}

#' @rdname is_null
#' @export
assert_is_null <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_null, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
