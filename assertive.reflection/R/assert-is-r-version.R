#' @rdname is_r
#' @export
assert_is_r_alpha <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_alpha, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_r_beta <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_beta, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_r_devel <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_devel, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_r_patched <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_patched, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_r_release_candidate <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_release_candidate, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_r_release <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_release, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_r_revised <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_r_revised, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_r_stable <- function(severity = getOption("assertive.severity", "stop"))
{                  
  .Deprecated("is_r_release")
  assert_engine(is_r_release, severity = severity)        
}

