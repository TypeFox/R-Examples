#' @rdname is_formula
#' @export
assert_is_formula <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_formula, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_formula
#' @export
assert_is_one_sided_formula <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_one_sided_formula, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_formula
#' @export
assert_is_two_sided_formula <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_two_sided_formula, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
