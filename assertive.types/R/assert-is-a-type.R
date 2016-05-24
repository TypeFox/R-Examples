#' @include imports.R
#' @rdname is_logical
#' @export
assert_is_a_bool <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{      
  assert_engine(
    is_a_bool, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

#' @rdname is_complex
#' @export
assert_is_a_complex <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                          
  assert_engine(
    is_a_complex, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

#' @rdname is_numeric
#' @export
assert_is_a_double <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                          
  assert_engine(
    is_a_double, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

#' @rdname is_numeric
#' @export
assert_is_a_number <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                          
  assert_engine(
    is_a_number, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

#' @rdname is_raw
#' @export
assert_is_a_raw <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                          
  assert_engine(
    is_a_raw, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  ) 
}

#' @rdname is_character
#' @export
assert_is_a_string <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_a_string, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

#' @rdname is_integer
#' @export
assert_is_an_integer <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    is_an_integer, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  ) 
}

#' @rdname is_inherited_from
#' @export
assert_is_inherited_from <- function(x, classes, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    is_inherited_from, 
    x, 
    classes = classes, 
    .xname = get_name_in_parent(x)
  )
}
