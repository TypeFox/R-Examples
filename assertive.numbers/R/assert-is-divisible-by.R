#' @rdname is_divisible_by
#' @export
assert_all_are_divisible_by <- function(x, n, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all divisible by %s (tol = %g).", 
    .xname, 
    toString(n, width = 20),
    tol
  )
  assert_engine(
    is_divisible_by, 
    x, 
    n = n, 
    tol = tol, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_divisible_by
#' @export
assert_any_are_divisible_by <- function(x, n, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are never divisible by %s (tol = %g).", 
    .xname, 
    toString(n, width = 20),
    tol
  )
  assert_engine(
    is_divisible_by, 
    x, 
    n = n, 
    tol = tol, 
    .xname = .xname,
    msg = msg,
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )   
}

#' @rdname is_divisible_by
#' @export
assert_all_are_even <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all even (tol = %g).", 
    .xname,
    tol
  )
  assert_engine(
    is_even, 
    x, 
    tol = tol,
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_divisible_by
#' @export
assert_any_are_even <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are never even (tol = %g).", 
    .xname,
    tol
  )
  assert_engine(
    is_even, 
    x, 
    tol = tol,
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  ) 
}

#' @rdname is_divisible_by
#' @export
assert_all_are_odd <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all odd (tol = %g).",
    .xname,
    tol
  )
  assert_engine(
    is_odd, 
    x, 
    tol = tol,
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_divisible_by
#' @export
assert_any_are_odd <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{  
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are never odd (tol = %g).", 
    .xname,
    tol
  )  
  assert_engine(
    is_odd, 
    x, 
    tol = tol,
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )
}

