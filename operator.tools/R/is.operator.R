#' Utilities for operators
#' 
#' These S4 Methods are utilies for working with operators. In R, operators are
#' functions with special syntax.
#' 
#' \code{is.operator} tests whether the object is one of the defined
#' \code{\link{operators}}.
#' 
#' \code{can.operator} tests whether the object can be coerced to an operator.
#' 
#' \code{as.operator} coerced the object to an operator.
#' 
#' Optionally, you can specify one of the that it tests for a specific type of
#' operator.  See details, below.
#' 
#' An operator is R function with special syntax.
#' 
#' ( See \code{??operator} for examples of each. )
#' 
#' \code{is.operator} tests whether the argumenst is an operator.
#' 
#' \code{as.operator} coerces \code{x} to a operator, otherwise fails.
#' 
#' \code{can.operator} test whether the object can be coerced to an operator.
#' 
#' All functions can accepts a \code{types} argument which is passed to
#' \code{link{operators}}.  By specifying one or more types, these functions
#' test using those \code{types} only.
#' 
#' New operators can be "registered" using \code{\link{setOperator}}.
#' 
#' @param x object to be tested or coerced. Can be \code{function} or
#' \code{name}.
#' @param \dots additional arguments passed to \code{\link{operators}}.
#' @return \code{is.operator} and \code{can.operator} return logical.
#' 
#' \code{as.operator} returns the argument coerced to the concommitant R
#' function.
#' @author Christopher Brown
#' @seealso \code{\link{operators}}, \code{\link{apropos}},
#' \code{\link{match.fun}}
#' @keywords utilities manip
#' @examples
#'  
#' 
#'  \dontrun{
#'    is.operator( `+` )
#'    is.operator( 'xyzzy' )
#'    is.operator( `+`, types="arithmetic" )
#'    is.operator( `+`, types="relational" )
#' 
#'    can.operator( `+` )
#'    can.operator( 'xyzzy' )
#'    can.operator( `+`, types="arithmetic" )
#'    can.operator( `+`, types="relational" )
#' 
#'    as.operator( `+` )
#'    as.operator( '+' )
#'    as.operator( as.name('+') )  
#'  }
#' 

#' @export 
is.operator <- function(x,...) 
  UseMethod( 'is.operator', x ) 

# DEFAULT: FALSE

#' @export
is.operator.default <- function(x, ... ) FALSE 

# NAME:
#' @export 
is.operator.name <- function(x, ... ) {
  x <- as.character(x) 

  if( length(list(...)) > 0 ) {
    x %in% operators(...)  
  } else {
    x %in% operators('REG') || x %in% operators('UNREG')
  }
   
}  

# FUNCTION: 
#' @export
is.operator.function <- function(x,...) 

  if( length(list(...)) > 0 ) {
    any( sapply(  
      operators(...), 
      function(op) identical( x, name2fun(op) )  
    ) ) 
  } else { 
    any( sapply( 
      operators("REG") ,
      function(op) identical( x, name2fun(op) )
    ) ) || 
    any( sapply( 
      operators("UNREG") ,
      function(op) identical( x, name2fun(op) )
    ) )  
  }      
