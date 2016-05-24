#' @rdname is_try_error
#' @export
assert_is_try_error <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_try_error, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_try_error
#' @export
assert_is_simple_error <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_simple_error, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_try_error
#' @export
assert_is_error <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_error, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_try_error
#' @export
assert_is_simple_warning <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_simple_warning, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_try_error
#' @export
assert_is_warning <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_warning, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_try_error
#' @export
assert_is_simple_message <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_simple_message, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_try_error
#' @export
assert_is_message <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_message, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_try_error
#' @export
assert_is_condition <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_warning, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

