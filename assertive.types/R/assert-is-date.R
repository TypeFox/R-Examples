#' @rdname is_date
#' @export
assert_is_date <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_date, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_date
#' @export
assert_is_posixct <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_posixct, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_date
#' @export
assert_is_posixlt <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_posixlt, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
