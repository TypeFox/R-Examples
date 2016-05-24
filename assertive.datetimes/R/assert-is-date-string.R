#' @rdname is_date_string
#' @export
assert_all_are_date_strings <- function(x, format = "%F %T", na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{    
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s contains strings that are not in the date format '%s'.", 
    get_name_in_parent(x),
    format, 
    .xname
  )
  assert_engine(
    is_date_string, 
    x, 
    format = format,
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_date_string
#' @export
assert_any_are_date_strings <- function(x, format = "%F %T", na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                   
  .xname <- get_name_in_parent(x)                                  
  msg <- gettextf(
    "%s contains no strings that are in the format '%s'.", 
    get_name_in_parent(x),
    format, 
    .xname
  )
  assert_engine(
    is_date_string, 
    x, 
    format = format,
    .xname = .xname,
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )      
}
