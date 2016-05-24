#' @rdname is_empty
#' @export
assert_has_elements <- function(x, n, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    has_elements, 
    x, 
    n = n,
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_empty
#' @export
assert_is_empty <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                             
  metric <- match.arg(metric)                             
  assert_engine(
    is_empty, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  ) 
}

#' @rdname is_empty
#' @export
assert_is_non_empty <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                            
  metric <- match.arg(metric)                                 
  assert_engine(
    is_non_empty, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  )  
}

#' @rdname is_empty
#' @export
assert_is_non_scalar <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                            
  metric <- match.arg(metric)                                 
  assert_engine(
    is_non_scalar, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  )   
}

##' @rdname is_empty
#' @export
assert_is_of_dimension <- function(x, n, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    is_of_dimension, 
    x, 
    n = n,
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_empty
#' @export
assert_is_of_length <- function(x, n, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    is_of_length, 
    x, 
    n = n,
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_empty
#' @export
assert_is_scalar <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                                        
  metric <- match.arg(metric)
  assert_engine(
    is_scalar, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
