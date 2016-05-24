#' Wrapper for writing S4 generics and methods
#' 
#' These are two wrappers around \code{setGeneric} and \code{setMethod}. A
#' relevant difference is that generics and methods are stored in the
#' environment in which  \code{\%g\%} and \code{\%m\%} are called and not in the
#' top-environment. Furthermore both functions have side effects in that they
#' will call \code{\link{globalVariables}} for the arguments and name of the
#' generic.
#' 
#' @param lhs see details
#' 
#' @param rhs the body as an expression
#' 
#' @details 
#' The Syntax for the left hand side:
#'   \cr\code{[<valueClass>:]<genericName>(<argList>)}
#'   \cr - \code{valueClass} optional, is the class of the return value (see
#'   \link{setGeneric})
#'   \cr - \code{genericName} the name of the generic function
#'   \cr - \code{argList} are \code{name = value} or \code{name ~ type}
#'   expressions. Name-Value expressions are just like in a function definition.
#'   Name-Type expressions are used to define the signature of a method (see
#'   \link{setMethod}). See \link{\%type\%} and the examples how to work with 
#'   them.
#' 
#' @examples 
#' # A new generic function and a method:
#' numeric : generic(x) %g% standardGeneric("generic") 
#' generic(x ~ numeric) %m% x
#' generic(1)
#' 
#' # Polymorphic methods in an object:
#' Object <- function() {
#'   numeric : generic(x) %g% standardGeneric("generic") 
#'   generic(x ~ numeric) %m% x
#'   retList("Object")
#' }
#' Object()$generic(1)
#' 
#' # Class Unions:
#' ## This generic allows for return values of type numeric or character:
#' 'numeric | character' : generic(x) %g% standardGeneric("generic")
#' 
#' ## This method also allows for numeric or character as argument:
#' generic(x ~ character | numeric) %m% x
#' generic(1)
#' generic("")
#' 
#' @export
#' @rdname S4generics
"%g%" <- function(lhs, rhs) {
  
  argList <- GenericExpressionTree(match.call(), parent.frame())
  
  # Fix for R CMD check. This can result in an error, if you use local generics
  # or methods.
  try(silent = TRUE, {
    globalVariables(c(
      argList$name, 
      names(formals(argList$def))), 
      topenv(argList$where)
    )
  })
  
  do.call("setGeneric", argList)
  invisible(getGeneric(argList$name, where = argList$where))
  
}

#' @export
#' @rdname S4generics
"%m%" <- function(lhs, rhs) {
  
  argList <- MethodExpressionTree(match.call(), parent.frame())
  
  # Fix for R CMD check. This can result in an error, if you use local generics
  # or methods.
  try(silent = TRUE, {
    globalVariables(
      names(formals(argList$definition)), 
      topenv(argList$where)
    )
  })
  
  do.call("setMethod", argList)
  invisible(getMethod(argList$f, argList$signature))
  
}
