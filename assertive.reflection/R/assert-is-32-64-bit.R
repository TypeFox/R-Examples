#' @include imports.R
#' 
#' @rdname is_windows
#' @export
assert_is_64_bit_os <- function(severity = getOption("assertive.severity", "stop"))
{
  .Deprecated("assert_is_64_bit")
  assert_engine(is_64_bit, severity = severity)
}

#' @rdname is_windows
#' @export
assert_is_32_bit <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(is_32_bit, severity = severity)
}

#' @rdname is_windows
#' @export
assert_is_64_bit <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(is_64_bit, severity = severity)
}
