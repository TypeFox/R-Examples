#' @rdname is_r
#' @export
assert_is_architect <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_architect, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_revo_r <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_revo_r, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_rstudio <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_rstudio, severity = severity)        
}

#' @rdname is_rstudio_current
#' @export
assert_is_rstudio_current <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_rstudio_current, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_rstudio_desktop <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_rstudio_desktop, severity = severity)        
}

#' @rdname is_r
#' @export
assert_is_rstudio_server <- function(severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(is_rstudio_server, severity = severity)        
}

