#' @rdname is_data.table
#' @export
assert_is_data.table <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_data.table, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

