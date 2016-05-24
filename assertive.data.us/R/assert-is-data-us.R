#' @include imports.R

#' @rdname is_us_social_security_number
#' @export
assert_all_are_us_social_security_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all US social security numbers.", get_name_in_parent(x))
  assert_engine(
    is_us_social_security_number, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_us_social_security_number
#' @export
assert_any_are_us_social_security_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not US social security numbers.", get_name_in_parent(x))
  assert_engine(
    is_us_social_security_number, 
    x, 
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )       
}

#' @rdname is_us_telephone_number
#' @export
assert_all_are_us_telephone_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all US telephone numbers.", get_name_in_parent(x))
  assert_engine(
    is_us_telephone_number, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_us_telephone_number
#' @export
assert_any_are_us_telephone_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not US telephone numbers.", get_name_in_parent(x))
  assert_engine(
    is_us_telephone_number, 
    x, 
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )       
}


#' @rdname is_us_zip_code
#' @export
assert_all_are_us_zip_codes <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all US zip codes.", get_name_in_parent(x))
  assert_engine(
    is_us_zip_code, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )         
}

#' @rdname is_us_zip_code
#' @export
assert_any_are_us_zip_codes <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not US zip codes.", get_name_in_parent(x))
  assert_engine(
    is_us_zip_code, 
    x, 
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )        
}
