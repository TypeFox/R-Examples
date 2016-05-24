#' @rdname is_r
#' @export
assert_is_r <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r, severity = severity)        
}

#' @rdname is_batch_mode
#' @export
assert_is_batch_mode <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(is_batch_mode, severity = severity)
}

#' @rdname is_batch_mode
#' @export
assert_is_interactive <- function(severity = getOption("assertive.severity", "stop"))
{
  assert_engine(is_interactive, severity = severity)
}

#' @rdname is_batch_mode
#' @export
assert_is_r_slave <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_slave, severity = severity)        
}

#' @rdname is_batch_mode
#' @export
assert_is_slave_r <- function(severity = getOption("assertive.severity", "stop"))
{                    
  .Deprecated("is_r_slave")                                     
  assert_engine(is_r_slave, severity = severity)        
}
