#' @rdname has_names
#' @export
assert_has_colnames <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                       
  assert_engine(
    has_colnames, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname has_names
#' @export
assert_has_dimnames <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                       
  assert_engine(
    has_dimnames, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname has_names
#' @export
assert_has_names <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                            
  assert_engine(
    has_names, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname has_names
#' @export
assert_has_rownames <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                             
  assert_engine(
    has_rownames, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
