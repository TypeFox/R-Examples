#' @rdname is_whole_number
#' @export
assert_all_numbers_are_whole_numbers <- function(x,
  tol = 100 * .Machine$double.eps, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                       
  .Deprecated("assert_all_are_whole_numbers")
  assert_all_are_whole_numbers(
    x, 
    tol = tol,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_whole_number
#' @export
assert_any_numbers_are_whole_numbers <- function(x, 
  tol = 100 * .Machine$double.eps, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                      
  .Deprecated("assert_any_are_whole_numbers")
  assert_any_are_whole_numbers(
    x, 
    tol = tol,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_whole_number
#' @export
assert_all_are_whole_numbers <- function(x,
  tol = 100 * .Machine$double.eps, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                               
  .xname <- get_name_in_parent(x)                   
  msg <- gettextf(
    "%s are not all whole numbers (tol = %g).", 
    .xname,
    tol
  )
  assert_engine(
    is_whole_number, 
    x, 
    tol = tol,
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore, 
    severity = severity
  ) 
}

#' @rdname is_whole_number
#' @export
assert_any_are_whole_numbers <- function(x, 
  tol = 100 * .Machine$double.eps, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{         
  .xname <- get_name_in_parent(x)                                      
  msg <- gettextf(
    "%s are never whole numbers (tol = %g).", 
   .xname,
    tol
  )
  assert_engine(
    is_whole_number, 
    x, 
    tol = tol,
    msg = msg,
    .xname = .xname,
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}
