#' @include imports.R

#' @rdname is_cas_number
#' @export
assert_all_are_cas_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all CAS numbers.", get_name_in_parent(x))
  assert_engine(
    is_cas_number, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )        
}

#' @rdname is_cas_number
#' @export
assert_any_are_cas_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not CAS numbers.", get_name_in_parent(x))
  assert_engine(
    is_cas_number, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  ) 
}

#' @rdname is_credit_card_number
#' @export
assert_all_are_credit_card_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all credit card numbers.", get_name_in_parent(x))
  assert_engine(
    is_credit_card_number, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  ) 
}

#' @rdname is_credit_card_number
#' @export
assert_any_are_credit_card_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not credit card numbers.", get_name_in_parent(x))
  assert_engine(
    is_credit_card_number, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )        
}

#' @rdname is_email_address
#' @export
assert_all_are_email_addresses <- function(x, method = c("simple", "rfc5322"), 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{                    
  method <- match.arg(method)
  msg <- gettextf("%s are not all email addresses.", get_name_in_parent(x))
  assert_engine(
    is_email_address, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_email_address
#' @export
assert_any_are_email_addresses <- function(x, method = c("simple", "rfc5322"), 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{                                 
  method <- match.arg(method)                    
  msg <- gettextf("%s are all not email addresses.", get_name_in_parent(x))
  assert_engine(
    is_email_address, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )       
}

#' @rdname is_hex_color
#' @export
assert_all_are_hex_colors <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                    
  msg <- gettextf("%s are not all hex colors.", get_name_in_parent(x))
  assert_engine(
    is_hex_color, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_hex_color
#' @export
assert_any_are_hex_colors <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                
  msg <- gettextf("%s are all not hex colors.", get_name_in_parent(x))
  assert_engine(
    is_hex_color, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )         
}

#' @rdname is_hex_color
#' @export
assert_all_are_hex_colours <- assert_all_are_hex_colors

#' @rdname is_hex_color
#' @export
assert_any_are_hex_colours <- assert_any_are_hex_colors

#' @rdname is_honorific
#' @export
assert_all_are_honorifics <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all honorifics.", get_name_in_parent(x))
  assert_engine(
    is_honorific, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )       
}

#' @rdname is_honorific
#' @export
assert_any_are_honorifics <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not honorifics.", get_name_in_parent(x))
  assert_engine(
    is_honorific, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_ip_address
#' @export
assert_all_are_ip_addresses <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all IP addresses.", get_name_in_parent(x))
  assert_engine(
    is_ip_address, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )       
}

#' @rdname is_ip_address
#' @export
assert_any_are_ip_addresses <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not IP addresses.", get_name_in_parent(x))
  assert_engine(
    is_ip_address, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )     
}

#' @rdname is_isbn_code
#' @export
assert_all_are_isbn_codes <- function(x, type = c("10", "13"), na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all ISBN codes.", get_name_in_parent(x))
  assert_engine(
    is_isbn_code, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )        
}

#' @rdname is_isbn_code
#' @export
assert_any_are_isbn_codes <- function(x, type = c("10", "13"), na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not ISBN codes.", get_name_in_parent(x))
  assert_engine(
    is_isbn_code, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )     
}
