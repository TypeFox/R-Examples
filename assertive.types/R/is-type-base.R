#' Is the input an array or matrix?
#'
#' Checks to see if the input is an array or matrix.
#'
#' @param x Input to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_array} and \code{is_matrix} wrap \code{is.array}, 
#' and \code{is.matrix} respectively, providing more information on
#' failure.  The \code{assert_*} functions return nothing but throw
#' an error if the corresponding \code{is_*} function returns
#' \code{FALSE}.
#' @examples
#' assert_is_array(array())
#' assert_is_array(matrix())
#' assert_is_matrix(matrix())
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_matrix(array()))
#' @importFrom assertive.base is2
#' @export
is_array <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "array", .xname)
}

#' @rdname is_language
#' @importFrom assertive.base is2
#' @export
is_call <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "call", .xname)
}

#' Is the input of type character?
#'
#' Checks to see if the input is of type character.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_character} wraps \code{is.character}, providing more 
#' information on failure. \code{is_a_string} returns \code{TRUE} if the 
#' input is character and scalar.
#' The \code{assert_*} functions return nothing but throw an error if the
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{is.character}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_character(letters)
#' assertive.base::dont_stop(assert_is_character(factor(letters)))
#' @importFrom assertive.base is2
#' @export
is_character <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "character", .xname)
}

#' Is the input complex?
#'
#' Checks to see if the input is complex.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_complex} wraps \code{is.complex}, providing more 
#' information on failure. \code{is_a_complex} returns \code{TRUE} if the 
#' input is complex and scalar.  The \code{assert_*} functions return
#' nothing but throw an error if the corresponding \code{is_*} function
#' returns \code{FALSE}.
#' @seealso \code{\link[base]{is.complex}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_complex(c(1i, 2i))
#' assert_is_a_complex(1i)
#' assert_is_a_complex(1 + 0i)
#' assert_is_a_complex(NA_complex_)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_complex(1:10))
#' assertive.base::dont_stop(assert_is_a_complex(c(1i, 2i)))
#' assertive.base::dont_stop(assert_is_a_complex(complex()))
#' @importFrom assertive.base is2
#' @export
is_complex <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "complex", .xname)
}       

#' Is the input is a data.frame?
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_data.frame} wraps \code{is.data.frame}, 
#' providing more information on failure.  \code{assert_is_data.frame} 
#' returns nothing but throws an error if \code{is_data.frame} 
#' returns \code{FALSE}.
#' @seealso \code{\link[base]{is.data.frame}}.
#' @examples
#' assert_is_data.frame(data.frame())
#' assert_is_data.frame(datasets::CO2)
#' @importFrom assertive.base is2
#' @export
is_data.frame <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "data.frame", .xname)
}

#' @rdname is_numeric
#' @importFrom assertive.base is2
#' @export
is_double <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "double", .xname)
}

#' Is the input an environment?
#'
#' Checks to see if the input is an environment.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_environment} wraps \code{is.environment}, providing more 
#' information on failure.  \code{assert_is_environment} returns nothing
#' but throws an error if \code{is_environment} returns \code{FALSE}.
#' @seealso \code{\link[base]{is.environment}}.
#' @examples
#' assert_is_environment(new.env())
#' assert_is_environment(globalenv())
#' assert_is_environment(baseenv())
#' @importFrom assertive.base is2
#' @export
is_environment <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "environment", .xname)
}

#' @rdname is_language
#' @export
is_expression <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "expression", .xname)
}

#' Is the input a factor?
#'
#' Checks to see if the input is an factor.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_factor} wraps \code{is.factor}, providing more 
#' information on failure.  \code{assert_is_factor} returns nothing
#' but throws an error if \code{is_factor} returns \code{FALSE}.
#' @seealso \code{\link[base]{is.factor}}.
#' @examples
#' assert_is_factor(factor(sample(letters, 10)))
#' @importFrom assertive.base is2
#' @export
is_factor <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "factor", .xname)
}

#' Is the input a function?
#'
#' Checks to see if the input is a function.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_function}, \code{is_primitive} and \code{is_stepfun}
#' wrap \code{is.function}, \code{is.primitive} and \code{is.stepfun} 
#' repsectively, providing more information on failure.  The 
#' \code{assert_*} functions return nothing but throw an error if the
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{is.function}}.
#' @examples
#' assert_is_function(sqrt)
#' assert_is_function(function(){})
#' @importFrom assertive.base is2
#' @export
is_function <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "function", .xname)
}

#' Is the input an integer?
#'
#' Checks to see if the input is an integer.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_integer} wraps \code{is.integer}, providing more 
#' information on failure. \code{is_an_integer} returns \code{TRUE} if the 
#' input is an integer and scalar.  The \code{assert_*} functions return 
#' nothing but throw an error if the corresponding \code{is_*} function
#' returns \code{FALSE}.
#' @seealso \code{\link[base]{is.integer}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_integer(1:10)
#' assert_is_an_integer(99L)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_integer(c(1, 2, 3)))
#' assertive.base::dont_stop(assert_is_an_integer(1:10))
#' assertive.base::dont_stop(assert_is_an_integer(integer()))
#' @importFrom assertive.base is2
#' @export
is_integer <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "integer", .xname)
}

#' Is the input a language object?
#'
#' Checks to see if the input is a language object.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_call}, \code{is_expression}, \code{is_language}, 
#' \code{is_name} and \code{is_symbol} wrap the corresponding \code{is.*}
#' functions, providing more information on failure.  The \code{assert_*}
#' functions return nothing but throw an error if the corresponding
#' \code{is_*} function returns \code{FALSE}.
#' @note \code{is_name} and \code{is_symbol} are different names for 
#' the same function.
#' @seealso \code{\link[base]{is.call}}, \code{\link[base]{is.expression}}
#' \code{\link[base]{is.language}} and \code{\link[base]{is.name}}.
#' @examples
#' a_call <- call("sin", "pi")
#' assert_is_call(a_call)
#' assert_is_language(a_call)
#' an_expression <- expression(sin(pi))
#' assert_is_expression(an_expression)
#' assert_is_language(an_expression)
#' a_name <- as.name("foo")
#' assert_is_name(a_name)
#' assert_is_language(a_name)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_language(function(){}))
#' @export
is_language <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.language(x)) 
  {
    return(
      false(
        gettext("%s is not a language object (name, call or expression)."), 
        .xname
      )
    )
  }
  TRUE
}

#' Is the input a list?
#'
#' Checks to see if the input is a list.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_list} wraps \code{is.list}, providing more 
#' information on failure.
#' @seealso \code{\link[base]{is.list}}.
#' @examples
#' assert_is_list(list(1,2,3))
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_list(1:10))
#' @importFrom assertive.base is2
#' @export
is_list <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "list", .xname)
}

#' Is the input logical?
#'
#' Checks to see if the input is logical.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_logical} wraps \code{is.logical}, providing more 
#' information on failure. \code{is_a_bool} returns \code{TRUE} if the 
#' input is logical and scalar.  The \code{assert_*} functions return
#' nothing but throw an error if the corresponding \code{is_*} function
#' returns \code{FALSE}.
#' @seealso \code{\link[base]{is.logical}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_logical(runif(10) > 0.5)
#' assert_is_a_bool(TRUE)
#' assert_is_a_bool(NA)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_logical(1))
#' assertive.base::dont_stop(assert_is_a_bool(c(TRUE, FALSE)))
#' assertive.base::dont_stop(assert_is_a_bool(logical()))
#' @importFrom assertive.base is2
#' @export
is_logical <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "logical", .xname)
}       

#' @rdname is_array
#' @importFrom assertive.base is2
#' @export
is_matrix <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "matrix", .xname)
}

#' @rdname is_language
#' @importFrom assertive.base is2
#' @export
is_name <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "name", .xname)
}

#' Is the input numeric?
#'
#' Checks to see if the input is numeric.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_numeric} wraps \code{is.numeric}, providing more 
#' information on failure. \code{is_a_number} returns \code{TRUE} if the 
#' input is numeric and scalar.  The \code{assert_*} functions return nothing
#' but throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @note \code{numeric} means either double or integer, inc this case.
#' @seealso \code{\link{is_integer}}, \code{\link[base]{is.numeric}} and 
#' \code{\link{is_scalar}}.
#' @examples
#' # "numeric" fns work on double or integers; 
#' assert_is_numeric(1:10)
#' 
#' # Here we check for length 1 as well as type
#' assert_is_a_number(pi)
#' assert_is_a_number(1L)
#' assert_is_a_number(NA_real_)
#' 
#' # "double" fns fail for integers.
#' assert_is_a_double(pi)
#' 
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_numeric(c(TRUE, FALSE)))
#' assertive.base::dont_stop(assert_is_a_number(1:10))
#' assertive.base::dont_stop(assert_is_a_number(numeric()))
#' assertive.base::dont_stop(assert_is_double(1:10))
#' @importFrom assertive.base is2
#' @export
is_numeric <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "numeric", .xname)
}

#' @rdname is_factor
#' @export
is_ordered <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_factor(x))) 
  {
    return(ok)
  }
  if(!is.ordered(x))
  {
    return(false(gettext("%s is not an ordered factor."), .xname))
  }
  TRUE
}

#' @rdname is_function
#' @export
is_primitive <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_function(x))) return(ok)
  if(!is.primitive(x))
  {
    return(false(gettext("%s is not a primitive function."), .xname))
  }
  TRUE
} 

#' Is the input a QR decomposition of a matrix?
#'
#' Checks to see if the input is a QR decomposition of a matrix.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_qr} wraps \code{is.qr}, providing more 
#' information on failure.  \code{assert_is_qr} returns nothing but
#' throws an error if \code{is_qr} returns \code{FALSE}.
#' @seealso \code{\link[base]{is.qr}}.
#' @examples
#' assert_is_qr(qr(matrix(rnorm(25), nrow = 5)))
#' @importFrom assertive.base is2
#' @export
is_qr <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "qr", .xname)
}

#' Is the input raw?
#'
#' Checks to see if the input is raw.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_raw} wraps \code{is.raw}, providing more 
#' information on failure. \code{is_a_raw} returns \code{TRUE} if the 
#' input is raw and scalar.  The \code{assert_*} functions return nothing but
#' throws an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[base]{is.raw}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_raw(as.raw(1:10))
#' assert_is_a_raw(as.raw(255))
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_raw(c(TRUE, FALSE)))
#' assertive.base::dont_stop(assert_is_a_raw(as.raw(1:10)))
#' assertive.base::dont_stop(assert_is_a_raw(raw()))
#' @importFrom assertive.base is2
#' @export
is_raw <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "raw", .xname)
}

#' Is the input an S4 object?
#'
#' Checks to see if the input is an S4 object.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_S4} wraps \code{isS4}, providing more information on 
#' failure.  \code{assert_is_S4} returns nothing but throws an error if
#' \code{is_S4} returns \code{FALSE}.
#' @seealso \code{\link[base]{isS4}}.
#' @examples
#' assert_is_s4(getClass("MethodDefinition"))
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_s4(1:10))
#' @export
is_s4 <- function(x, .xname = get_name_in_parent(x))
{
  if(!isS4(x))
  {
    return(false(gettext("%s is not an S4 object."), .xname))
  }
  TRUE
} 

#' @rdname is_s4
#' @export
is_S4 <- function(x, .xname = get_name_in_parent(x))
{
  .Deprecated("is_s4")
  is_s4(x, .xname)
}

#' @rdname is_language
#' @export
is_symbol <- is_name

#' Is the input a table?
#'
#' Checks to see if the input is a table.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_table} wraps \code{is.table}, providing more 
#' information on failure.  \code{assert_is_table} returns nothing but
#' throws an error if \code{is_table} returns \code{FALSE}.
#' @seealso \code{\link[base]{is.table}}.
#' @examples
#' assert_is_table(table(sample(letters, 100, replace = TRUE)))
#' @importFrom assertive.base is2
#' @export
is_table <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "table", .xname)
}
