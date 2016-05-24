#' @rdname is_unsorted
#' @export
assert_is_unsorted <- function(x, na.rm = FALSE, strictly = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_unsorted, 
    x,
    na.rm = na.rm,
    strictly = strictly, 
    .xname = get_name_in_parent(x),
    severity = severity
  )       
}
