#' @rdname has_slot
#' @export
assert_has_slot <- function(x, severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    has_slot,
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
