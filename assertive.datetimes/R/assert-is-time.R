#' @include imports.R

#' @rdname is_in_past
#' @export
assert_all_are_after <- function(x, y, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{              
  .xname = get_name_in_parent(x)  
  .yname = get_name_in_parent(y)                                       
  msg <- gettextf("%s are not all after %s.", .xname, .yname)
  assert_engine(
    is_after,
    x, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  ) 
}

#' @rdname is_in_past
#' @export
assert_any_are_after <- function(x, y, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                
  .xname = get_name_in_parent(x)  
  .yname = get_name_in_parent(y)                                     
  msg <- gettextf("%s are all before %s.", .xname, .yname)
  assert_engine(
    is_after,
    x, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )   
}
 
#' @rdname is_in_past
#' @export
assert_all_are_before <- function(x, y, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{            
  .xname = get_name_in_parent(x)  
  .yname = get_name_in_parent(y)                                           
  msg <- gettextf("%s are not all before %s.", .xname, .yname)
  assert_engine(
    is_before,
    x, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  ) 
}

#' @rdname is_in_past
#' @export
assert_any_are_before <- function(x, y, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{            
  .xname = get_name_in_parent(x)  
  .yname = get_name_in_parent(y)                                       
  msg <- gettextf("%s are all after %s.", .xname, , .yname)
  assert_engine(
    is_before,
    x, 
    .xname = .xname,
    .yname = .yname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )   
}
 
#' @rdname is_in_past
#' @export
assert_all_are_in_future <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                            
  .xname = get_name_in_parent(x)                         
  msg <- gettextf("%s are not all in the future.", .xname)
  assert_engine(
    is_in_future,
    x, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  ) 
}

#' @rdname is_in_past
#' @export
assert_any_are_in_future <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                        
  .xname = get_name_in_parent(x)                              
  msg <- gettextf("%s are all in the past.", .xname)
  assert_engine(
    is_in_future, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )   
}
 
#' @rdname is_in_past
#' @export
assert_all_are_in_past <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                 
  .xname = get_name_in_parent(x)                                     
  msg <- gettextf("%s are not all in the past.", .xname)
  assert_engine(
    is_in_past, 
    x, 
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )    
}

#' @rdname is_in_past
#' @export
assert_any_are_in_past <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                 
  .xname = get_name_in_parent(x)                                     
  msg <- gettextf("%s are all in the future.", .xname)
  assert_engine(
    is_in_past, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )  
}
