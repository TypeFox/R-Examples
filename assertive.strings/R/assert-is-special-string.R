#' @rdname is_numeric_string
#' @export
assert_all_are_numeric_strings <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{        
  .xname <- get_name_in_parent(x)                                             
  msg <- gettextf(
    "%s is not a character vector of numbers.", 
    .xname
  )
  assert_engine(
    is_numeric_string, 
    x, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_numeric_string
#' @export
assert_any_are_numeric_strings <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{        
  .xname <- get_name_in_parent(x)                                             
  msg <- gettextf(
    "%s is not a character vector of numbers.", 
    .xname
  )
  assert_engine(
    is_numeric_string, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )         
}

#' @rdname is_numeric_string
#' @export
assert_all_are_logical_strings <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{        
  .xname <- get_name_in_parent(x)                                             
  msg <- gettextf(
    "%s is not a character vector of logical values.", 
    .xname
  )
  assert_engine(
    is_logical_string, 
    x, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_numeric_string
#' @export
assert_any_are_logical_strings <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{        
  .xname <- get_name_in_parent(x)                                             
  msg <- gettextf(
    "%s is not a character vector of logical values.", 
    .xname
  )
  assert_engine(
    is_logical_string, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )         
}

#' @rdname is_single_character
#' @export
assert_all_are_single_characters <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{         
  .xname <- get_name_in_parent(x)                                            
  msg <- gettextf("%s are not all single characters.", .xname)
  assert_engine(
    is_single_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_single_character
#' @export
assert_any_are_single_characters <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{            
  .xname <- get_name_in_parent(x)                                         
  msg <- gettextf("%s are all not single characters.", .xname)
  assert_engine(
    is_single_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )
}
