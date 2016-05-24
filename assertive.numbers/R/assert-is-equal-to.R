#' @rdname is_equal_to
#' @export
assert_all_are_equal_to <- function(x, y, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are not all equal to %s (tol = %g).", 
    .xname,
    .yname,
    tol
  )
  assert_engine(
    is_equal_to, 
    x, 
    y = y, 
    tol = tol, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_any_are_equal_to <- function(x, y, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are never equal to %s (tol = %g).", 
    .xname,
    .yname,
    tol
  )
  assert_engine(
    is_equal_to, 
    x, 
    y = y, 
    tol = tol, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_all_are_not_equal_to <- function(x, y, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are sometimes equal to %s (tol = %g).", 
    .xname,
    .yname,
    tol
  )
  assert_engine(
    is_not_equal_to, 
    x, 
    y = y, 
    tol = tol, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_any_are_not_equal_to <- function(x, y, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are all equal to %s (tol = %g).", 
    .xname,
    .yname,
    tol
  )
  assert_engine(
    is_not_equal_to, 
    x, 
    y = y, 
    tol = tol, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_all_are_greater_than <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are not all greater than %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_greater_than, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_any_are_greater_than <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are never greater than %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_greater_than, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_all_are_greater_than_or_equal_to <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are not all greater than or equal to %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_greater_than_or_equal_to, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_any_are_greater_than_or_equal_to <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are never greater than or equal to %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_greater_than_or_equal_to, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_all_are_less_than <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are not all less than %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_less_than, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_any_are_less_than <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are never less than %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_less_than, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_all_are_less_than_or_equal_to <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are not all less than or equal to %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_less_than_or_equal_to, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_equal_to
#' @export
assert_any_are_less_than_or_equal_to <- function(x, y, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  .yname <- get_name_in_parent(y)
  msg <- gettextf(
    "%s are never less than or equal to %s.", 
    .xname,
    .yname
  )
  assert_engine(
    is_less_than_or_equal_to, 
    x, 
    y = y, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )  
}
