#' @rdname is_finite
#' @export
assert_all_are_finite <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                  
  .xname <- get_name_in_parent(x)                                   
  msg <- gettextf("%s are not all finite.", .xname)
  assert_engine(
    is_finite, 
    x,
    .xname = .xname,
    msg = msg,
    severity = severity
  )        
}

#' @rdname is_finite
#' @export
assert_any_are_finite <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                   
  .xname <- get_name_in_parent(x)                                      
  msg <- gettextf("%s are never finite.", .xname)
  assert_engine(
    is_finite, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )        
}

#' @rdname is_finite
#' @export
assert_all_are_infinite <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)                                   
  msg <- gettextf("%s are not all infinite.", .xname)
  assert_engine(
  is_infinite, 
    x,
    .xname = .xname,
    msg = msg, 
    severity = severity
  )    
}

#' @rdname is_finite
#' @export
assert_any_are_infinite <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                           
  .xname <- get_name_in_parent(x)                          
  msg <- gettextf("%s are never infinite.", .xname)
  assert_engine(
    is_infinite, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )     
}

#' @rdname is_nan
#' @export
assert_all_are_nan <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                         
  .xname <- get_name_in_parent(x)                                       
  msg <- gettextf("%s are not all NaN.", .xname)
  assert_engine(
    is_nan, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_nan
#' @export
assert_any_are_nan <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                      
  .xname <- get_name_in_parent(x)                                          
  msg <- gettextf("%s are never NaN.", .xname)
  assert_engine(
    is_nan, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_finite
#' @export
assert_all_are_negative_infinity <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                       
  .xname <- get_name_in_parent(x)                              
  msg <- gettextf("%s are not all negative infinity.", .xname)
  assert_engine(
    is_negative_infinity, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )       
}

#' @rdname is_finite
#' @export
assert_any_are_negative_infinity <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                        
  .xname <- get_name_in_parent(x)                             
  msg <- gettextf("%s are never negative infinity.", .xname)
  assert_engine(
    is_negative_infinity, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )      
}

#' @rdname is_nan
#' @export
assert_all_are_not_nan <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                       
  .xname <- get_name_in_parent(x)                               
  msg <- gettextf("%s contains NaNs.", .xname)
  assert_engine(
    is_not_nan, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_nan
#' @export
assert_any_are_not_nan <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                    
  .xname <- get_name_in_parent(x)                                  
  msg <- gettextf("%s are all NaN.", .xname)
  assert_engine(
    is_not_nan, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_finite
#' @export
assert_all_are_positive_infinity <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                    
  .xname <- get_name_in_parent(x)                                 
  msg <- gettextf("%s are not all positive infinity.", .xname)
  assert_engine(
    is_positive_infinity, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )     
}

#' @rdname is_finite
#' @export
assert_any_are_positive_infinity <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                   
  .xname <- get_name_in_parent(x)                                  
  msg <- gettextf("%s are never positive infinity.",.xname)
  assert_engine(
    is_positive_infinity, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )     
}
