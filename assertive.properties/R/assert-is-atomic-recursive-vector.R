#' @rdname is_atomic
#' @export
assert_is_atomic <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                    
  assert_engine(
    is_atomic, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_atomic
#' @export
assert_is_nested <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_nested, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_atomic
#' @export
assert_is_non_nested <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_non_nested, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_atomic
#' @export
assert_is_recursive <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_recursive, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_atomic
#' @export
assert_is_vector <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                    
  assert_engine(
    is_vector, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
