#' @include imports.R

#' @rdname is_in_range
#' @export
assert_all_are_in_closed_range <- function(x, lower = -Inf, upper = Inf, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)                                          
  msg <- gettextf(
    "%s are not all in the range %s.", 
    .xname,
    make_range_string(lower, upper, FALSE, FALSE)
  )
  assert_engine(
    is_in_closed_range, 
    x,  
    lower = lower, 
    upper = upper,
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_any_are_in_closed_range <- function(x, lower = -Inf, upper = Inf, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{                      
  .xname <- get_name_in_parent(x)                                
  msg <- gettextf(
    "%s are all out of the range %s.", 
    .xname,
    make_range_string(lower, upper, FALSE, FALSE)
  )
  assert_engine(
    is_in_closed_range, 
    x, 
    lower = lower, 
    upper = upper, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_all_are_in_left_open_range <- function(x, lower = -Inf, upper = Inf, 
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{                   
  .xname <- get_name_in_parent(x)                                   
  msg <- gettextf(
    "%s are not all in the range %s.", 
    .xname,
    make_range_string(lower, upper, TRUE, FALSE)
  )
  assert_engine(
    is_in_left_open_range, 
    x, 
    lower = lower, 
    upper = upper, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_any_are_in_left_open_range <- function(x, lower = -Inf, upper = Inf,  
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{                   
  .xname <- get_name_in_parent(x)                                   
  msg <- gettextf(
    "%s are all out of the range %s.", 
    .xname,
    make_range_string(lower, upper, TRUE, FALSE)
  )
  assert_engine(
    is_in_left_open_range, 
    x, 
    lower = lower, 
    upper = upper,
    .xname = .xname,
    msg   = msg, 
    what  = "any",
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_all_are_in_open_range <- function(x, lower = -Inf, upper = Inf,  
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{              
  .xname <- get_name_in_parent(x)                                        
  msg <- gettextf(
    "%s are not all in the range %s.", 
    .xname,
    make_range_string(lower, upper, TRUE, TRUE)
  )
  assert_engine(
    is_in_open_range, 
    x, 
    lower = lower, 
    upper = upper, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_any_are_in_open_range <- function(x, lower = -Inf, upper = Inf,  
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{              
  .xname <- get_name_in_parent(x)                                        
  msg <-   msg <- gettextf(
    "%s are all out of the range %s.", 
    .xname,
    make_range_string(lower, upper, TRUE, TRUE)
  )
  assert_engine(
    is_in_open_range, 
    x,
    lower = lower, 
    upper = upper, 
    .xname = .xname,
    msg   = msg, 
    what  = "any",
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_all_are_in_range <- function(x, lower = -Inf, upper = Inf, 
  lower_is_strict = FALSE, upper_is_strict = FALSE,  
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{               
  .xname <- get_name_in_parent(x)                                       
  msg <- gettextf(
    "%s are not all in the range %s.", 
    .xname,
    make_range_string(lower, upper, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_in_range, 
    x, 
    lower = lower, 
    upper = upper, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  ) 
}

#' @rdname is_in_range
#' @export
assert_any_are_in_range <- function(x, lower = -Inf, upper = Inf, 
  lower_is_strict = FALSE, upper_is_strict = FALSE,  
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                         
  msg <- gettextf(
    "%s are all out of the range %s.", 
    .xname,
    make_range_string(lower, upper, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_in_range, 
    x, 
    lower = lower, 
    upper = upper, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_all_are_in_right_open_range <- function(x, lower = -Inf, upper = Inf,  
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                         
  msg <- gettextf(
    "%s are not all in the range %s.", 
    .xname,
    make_range_string(lower, upper, FALSE, TRUE)
  )
  assert_engine(
    is_in_right_open_range, 
    x, 
    .xname = .xname,
    msg = msg,
    lower = lower, 
    upper = upper,
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_any_are_in_right_open_range <- function(x, lower = -Inf, upper = Inf,  
  na_ignore = FALSE, severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                         
  msg <- gettextf(
    "%s are all out of the range %s.", 
    .xname,
    make_range_string(lower, upper, FALSE, TRUE)
  )
  assert_engine(
    is_in_right_open_range, 
    x, 
    .xname = .xname,
    msg = msg,
    lower = lower, 
    upper = upper, 
    .xname = .xname,
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )  
}

#' @rdname is_in_range
#' @export
assert_all_are_negative <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{            
  .xname <- get_name_in_parent(x)                                                     
  msg <- gettextf("%s are not all negative.", .xname)
  assert_engine(
    is_negative, 
    x, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_any_are_negative <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{           
  .xname <- get_name_in_parent(x)                                              
  msg <- gettextf("%s are never negative.", .xname)
  assert_engine(
    is_negative, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_all_are_non_negative <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{           
  .xname <- get_name_in_parent(x)                                             
  msg <- gettextf("%s are not all non-negative.", .xname)
  assert_engine(
    is_non_negative, 
    x, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_any_are_non_negative <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{           
  .xname <- get_name_in_parent(x)                                            
  msg <- gettextf("%s are all negative.", .xname)
  assert_engine(
    is_non_negative, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_all_are_non_positive <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                           
  msg <- gettextf("%s contains positive values.", .xname)
  assert_engine(
    is_non_positive, 
    x, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_any_are_non_positive <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                          
  msg <- gettextf("%s are all positive.", .xname)
  assert_engine(
    is_non_positive, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_all_are_percentages <- function(x, lower_is_strict = FALSE, 
  upper_is_strict = FALSE, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{             
  .xname <- get_name_in_parent(x)                                           
  msg <- gettextf(
    "%s are not all in the range %s.", 
    .xname,
    make_range_string(0, 100, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_percentage, 
    x, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_any_are_percentages <- function(x, lower_is_strict = FALSE, 
  upper_is_strict = FALSE, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{           
  .xname <- get_name_in_parent(x)                                             
  msg <- gettextf(
    "%s are all out of the range %s.", 
    .xname,
    make_range_string(0, 100, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_percentage, 
    x, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}  

#' @rdname is_in_range
#' @export
assert_all_are_positive <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{           
  .xname <- get_name_in_parent(x)                                             
  msg <- gettextf("%s contains non-positive values.", .xname)
  assert_engine(
    is_positive, 
    x, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_any_are_positive <- function(x, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{           
  .xname <- get_name_in_parent(x)                                            
  msg <- gettextf("%s are all non-positive.", .xname)
  assert_engine(
    is_positive, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_all_are_proportions <- function(x, lower_is_strict = FALSE, 
  upper_is_strict = FALSE, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{        
  .xname <- get_name_in_parent(x)                                                
  msg <- gettextf(
    "%s are not all in the range %s.", 
    .xname,
    make_range_string(0, 1, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_proportion, 
    x, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict, 
    .xname = .xname,
    msg = msg,
    na_ignore = na_ignore, 
    severity = severity
  )
}

#' @rdname is_in_range
#' @export
assert_any_are_proportions <- function(x, lower_is_strict = FALSE, 
  upper_is_strict = FALSE, na_ignore = FALSE, 
severity = getOption("assertive.severity", "stop"))
{       
  .xname <- get_name_in_parent(x)                                                 
  msg <- gettextf(
    "%s are all out of the range %s.", 
    .xname,
    make_range_string(0, 1, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_proportion, 
    x, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore, 
    severity = severity
  )
}  

#' Make a range string
#' 
#' Makes a range string using mathematical notation.
#' @param lower A number for the lower bound.
#' @param upper A number for the upper bound.
#' @param lower_is_strict Should the lower bound be included?
#' @param upper_is_strict Should the upper bound be included?
#' @return A string denoting the range.
#' @note Not vectorized across the \code{lower_is_strict} and 
#' \code{upper_is_strict} args for speed.
#' @examples 
#' \donttest{
#' make_range_string(-1.2345, 6.7890, TRUE, FALSE)
#' }
#' @noRd
make_range_string <- function(lower, upper, lower_is_strict, upper_is_strict)
{
  left <- if(lower_is_strict) "(" else "["
  right <- if(upper_is_strict) ")" else "]"
  paste0(left, lower, ",", upper, right)  
}
