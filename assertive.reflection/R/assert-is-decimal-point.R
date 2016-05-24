#' @rdname is_xxx_for_decimal_point
#' @export
assert_is_comma_for_decimal_point <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(is_comma_for_decimal_point, severity = severity)
}

#' @rdname is_xxx_for_decimal_point
#' @export
assert_is_period_for_decimal_point <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(is_period_for_decimal_point, severity = severity)
}
