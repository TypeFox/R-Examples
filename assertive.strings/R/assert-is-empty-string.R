#' @rdname is_empty_character
#' @export
assert_is_an_empty_string <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    is_an_empty_string, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )    
}

#' @rdname is_empty_character
#' @export
assert_is_a_non_empty_string <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  assert_engine(
    is_a_non_empty_string, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )    
}  

#' @rdname is_empty_character
#' @export
assert_is_a_missing_or_empty_string <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                  
  assert_engine(
    is_a_missing_or_empty_string, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  ) 
}

#' @rdname is_empty_character
#' @export
assert_is_a_non_missing_nor_empty_string <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  assert_engine(
    is_a_non_missing_nor_empty_string, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}  
