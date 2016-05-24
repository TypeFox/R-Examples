#' @include imports.R
#' 
#' @rdname has_arg
#' @export
assert_has_arg <- function(x, fn = sys.function(sys.parent()), 
  severity = getOption("assertive.severity", "stop"))
{
  xname <- get_name_in_parent(x)
  msg <- gettextf(
    "The function %s does not have the argument %s.", 
    get_name_in_parent(fn), 
    xname
  )
  assert_engine(
    has_arg_,
    xname,
    fn = fn,
    msg = msg,
    severity = severity
  )  
}

#' @rdname is_binding_locked
#' @export
assert_is_binding_locked <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{      
  msg <- gettextf("%s does not have a locked binding.", get_name_in_parent(x))
  assert_engine(
    is_binding_locked,
    x,
    msg = msg,
    severity = severity
  )    
}

#' @rdname is_debugged
#' @export
assert_is_debugged <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{             
  msg <- gettextf("%s is not a function being debugged.", get_name_in_parent(x))
  assert_engine(
    is_debugged, 
    x, 
    msg = msg,
    severity = severity
  )       
}

#' @rdname is_existing
#' @export
assert_all_are_existing <- function(
  x, 
  envir = parent.frame(),
  inherits = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{    
  msg <- gettextf("%s do not all exist.", get_name_in_parent(x))
  assert_engine(
    is_existing, 
    x,
    envir = envir,
    inherits = inherits, 
    msg = msg,
    severity = severity
  )       
}

#' @rdname is_existing
#' @export
assert_any_are_existing <- function(
  x, 
  envir = parent.frame(), 
  inherits = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{    
  msg <- gettextf("%s all do not exist.", get_name_in_parent(x))
  assert_engine(
    is_existing, 
    x, 
    envir = envir,
    inherits = inherits,
    msg = msg,
    what = "any",
    severity = severity
  )       
}

# ' @rdname is_generic
# ' @export
# assert_is_generic <- function(x)
# {                                                     
#   msg <- gettextf("%s is not a generic function.", get_name_in_parent(x))
#   assert_engine(is_generic, x, msg = msg)        
# }

#' @rdname is_if_condition
#' @export
assert_is_if_condition <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    is_if_condition, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_loaded
#' @export
assert_is_loaded <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{          
  msg <- gettextf("%s is not loaded.", get_name_in_parent(x))                                               
  assert_engine(
    is_loaded, 
    x, 
    msg = msg,
    severity = severity
  )  
}

#' @rdname is_valid_r_code
#' @export
assert_is_valid_r_code <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                      
  msg <- gettextf("%s is not valid R code.", get_name_in_parent(x))
  assert_engine(
    is_valid_r_code, 
    x, 
    msg = msg,
    severity = severity
  )
}

#' @rdname is_valid_variable_name
#' @export
assert_all_are_valid_variable_names <- function(x, allow_reserved = TRUE, 
  allow_duplicates, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{   
  if(!missing(allow_duplicates))
  {
    .Deprecated(
      msg = "The 'allow_duplicates' argument is deprecated and will be ignored."
    )
  }
  msg <- gettextf("%s are not all valid variable names.", get_name_in_parent(x))
  assert_engine(
    is_valid_variable_name, 
    x,
    allow_reserved = allow_reserved, 
    msg = msg,
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_valid_variable_name
#' @export
assert_any_are_valid_variable_names <- function(x, allow_reserved = TRUE, 
  allow_duplicates, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{          
  if(!missing(allow_duplicates))
  {
    .Deprecated(
      msg = "The 'allow_duplicates' argument is deprecated and will be ignored."
    )
  }
  msg <- gettextf("%s are all invalid variable names.", get_name_in_parent(x))
  assert_engine(
    is_valid_variable_name, 
    x,
    allow_reserved = allow_reserved, 
    msg = msg,
    what = "any",
    na_ignore = na_ignore,
    severity = severity
  )
}
