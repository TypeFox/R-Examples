#' can.operator
#' 
#' tests whether an object can be coerced to a operator, optionally an operator
#' of 'types'.
#' 
#' 
#' @param x object; to test
#' @param ... additional arguments 
#' 
#' \code{can.operator} test whether an object can be coerced to an operator.
#' Methods exist for name, function amd character classes
#' 
#' @return logical

#' @export
can.operator <- function(x, ... )
  UseMethod( 'can.operator', x ) 

#' @export
can.operator.default <- function(x,...) FALSE

#' @export
can.operator.name <- function(x,...) 
  deparse(x) %in% operators(...)

#' @export
can.operator.function <- function(x,...) {
  deparse( x ) %in% 
  lapply( operators(...), function(op) deparse( eval( as.symbol(op) ) ) ) 
}  

#' @export
can.operator.character <- function(x,...) 
   x %in% operators(...)
