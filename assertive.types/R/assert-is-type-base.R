#' Does x belong to these classes?
#' 
#' Checks to see if x belongs to any of the classes in classes.
#' @param x Input to check.
#' @param classes As for \code{class}. 
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The functions return nothing but throw an error if
#' \code{x} does not have any/all of the class \code{classes}.
#' @seealso \code{\link[assertive.base]{is2}}
#' @examples 
#' assert_is_all_of(1:10, c("integer", "numeric"))
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_any_of(1:10, c("list", "data.frame")))
#' @export
assert_is_all_of <- function(x, classes, 
  severity = getOption("assertive.severity", "stop"))
{  
  msg <- gettextf(
    "%s is not in all of the classes %s.", 
    get_name_in_parent(x), 
    toString(sQuote(classes))
  )
  assert_engine(
    is2, 
    x, 
    class = classes, 
    msg = msg, 
    severity = severity
  )
}

#' @rdname assert_is_all_of  
#' @export
assert_is_any_of <- function(x, classes, 
  severity = getOption("assertive.severity", "stop"))
{  
  msg <- gettextf(
    "%s is not in any of the classes %s.", 
    get_name_in_parent(x), 
    toString(sQuote(classes))
  )
  assert_engine(
    is2, 
    x, 
    class = classes, 
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_array
#' @export
assert_is_array <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_array, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_language
#' @export
assert_is_call <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_call, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_character
#' @export
assert_is_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_character, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_complex
#' @export
assert_is_complex <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_complex, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_data.frame
#' @export
assert_is_data.frame <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_data.frame, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_numeric
#' @export
assert_is_double <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_double, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_environment
#' @export
assert_is_environment <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_environment, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_language
#' @export
assert_is_expression <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_expression, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_factor
#' @export
assert_is_factor <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_factor, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_function
#' @export
assert_is_function <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_function, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_integer
#' @export
assert_is_integer <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_integer, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_language
#' @export
assert_is_language <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_language, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_list
#' @export
assert_is_list <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_list, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_logical
#' @export
assert_is_logical <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_logical, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_array
#' @export
assert_is_matrix <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_matrix, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_language
#' @export
assert_is_name <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_name, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_numeric
#' @export
assert_is_numeric <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_numeric, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_factor
#' @export
assert_is_ordered <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_ordered, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_function
#' @export
assert_is_primitive <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_primitive, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_qr
#' @export
assert_is_qr <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_qr, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_raw
#' @export
assert_is_raw <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_raw, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_s4
#' @export
assert_is_S4 <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  .Deprecated("assert_is_s4")
  assert_is_s4(x)
}

#' @rdname is_s4
#' @export
assert_is_s4 <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_s4, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_language
#' @export
assert_is_symbol <- assert_is_name

#' @rdname is_table
#' @export
assert_is_table <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_table, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
