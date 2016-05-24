#' @rdname is_real
#' @export
assert_all_are_imaginary <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all imaginary (tol = %g).", 
    .xname, 
    tol
  )
  assert_engine(
    is_imaginary, 
    x, 
    tol = tol, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore, 
    severity = severity
  )        
}

#' @rdname is_real
#' @export
assert_any_are_imaginary <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{               
  .xname <- get_name_in_parent(x)                                          
  msg <- gettextf(
    "%s are never imaginary (tol = %g).", 
    .xname, 
    tol
  )
  assert_engine(
    is_imaginary, 
    x, 
    tol = tol, 
    .xname = .xname,
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore, 
    severity = severity
  )        
}

#' @rdname is_real
#' @export
assert_all_are_real <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                 
  msg <- gettextf(
    "%s are not all real (tol = %g).", 
    .xname, 
    tol
  )
  assert_engine(
    is_real, 
    x, 
    tol = tol, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_real
#' @export
assert_any_are_real <- function(x, tol = 100 * .Machine$double.eps, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                        
  msg <- gettextf("%s are never real (tol = %g).", .xname, tol)
  assert_engine(
    is_real, 
    x, 
    tol = tol,
    .xname = .xname,
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore, 
    severity = severity
  )
}
