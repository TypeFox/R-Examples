#' Types
#' 
#' This function can be used to define new S4-classes which are called Type.
#' They have an initialize method and in the introduced syntax init-method and
#' S4-class definition build a unit, hence a type. This simply captures a
#' typical \code{setClass} then \code{setMethod("initialize", ...)} pattern
#' where often some redundancy is introduced. The function has side effects due
#' to calling \code{setClass}, \code{setMethod} and assigning the constructor
#' function to the types name.
#' 
#' @param lhs an expression of the form:
#'   \cr\code{[<parent-name>:]<type-name>([<slots>])}
#'   \cr - <parent-name> optional, the name of the S4-class/type to inherit from, 
#' seperated by \code{:}
#'   \cr - <type-name> the name for the new type and constructor function.
#'   \cr - <slots> optional, \code{name = value} or \code{name ~ type}
#' expressions. Name-Value expressions are used to construct a prototype. From
#' the prototype the class of the slot will be inferred. They are also the
#' defaults in the type constructor. Name-Type expressions define the classes of
#' the slots. If no value (or type) is supplied, \code{ANY} is assumed.
#' @param rhs the body of the initialize method as expression. It will be called
#'   with \code{.Object} and \code{...} as arguments. \code{.Object} should be
#'   the return value. With \code{.Object} there is an instance of the type on
#'   which assertions can be formulated. Prior to the body (rhs) \code{.Object
#'   <- callNextMethod()} will be evaluated which enables proper initialization
#'   of your type and its inherited fields. See \link[methods]{initialize} for
#'   details.
#' 
#' @details 
#' \code{Name-Type} expressions are also used in \link{\%m\%}. Besides this you
#' can formulate type unions in type expressions or the inheritance structure.
#' This has a side effect in that \link{setClassUnion} is called. Whenever you
#' write a type you can replace the name by an expression of the form:
#' \code{type1 | type2}. Outside the slots or argument list of a method these
#' expressions have to be quoted. In this example the following expression is
#' evaluated for you: \code{setClassUnion("type1ORtype2", c("type1", "type2"))}.
#' 
#' 
#' @examples 
#' # This will create an S4-class named 'Test' with two slots; x = "numeric"
#' # and y = "list"; prototype: list(x = 1, y = list()); and an initialize
#' # method where some checks are performed.
#' 
#' Test(x = 1, y = list()) %type% {
#'   stopifnot(.Object@@x > 0)
#'   .Object
#' }
#' 
#' # This will create an S4-class named 'Numeric' with a slot and some tests.
#' 
#' numeric : Numeric(metaInfo = character()) %type% {
#'   stopifnot(length(.Object) > 0)
#'   stopifnot(all(.Object > 0))
#'   .Object
#' }
#' 
#' # This will create an S4-class with slots, where the constructor function has
#' # no defaults. All slots will allow for ANY type.
#' 
#' Anything(x, y ~ ANY, z = NULL) %type% .Object
#' \dontrun{
#'   Anything() # error because x and y are missing
#' }
#' 
#' # Type Unions:
#' 'character | numeric' : Either(either ~ character | numeric) %type% .Object
#' Either("", 1)
#' 
#' @export
"%type%" <- function(lhs, rhs) {
  
  .exprTree <- ExpressionTree(match.call(), parent.frame())
  .classExprTree <- ClassExpressionTree(.exprTree, parent.frame())
  .initExprTree <- InitMethodExpressionTree(.exprTree, parent.frame())
  .constExprTree <- ConstExpressionTree(.exprTree, parent.frame())
  
  do.call(setClass, .classExprTree)
  do.call(setMethod, .initExprTree)
  const <- do.call(makeFunDef, .constExprTree)
  assign(.exprTree$names[1], const, envir = parent.frame())
  
  globalVariables(c(
    .exprTree$names[1], 
    names(.classExprTree$slots)), 
    package = topenv(parent.frame())
  )
  
  invisible(getClass(.exprTree$names[1], where = parent.frame()))
  
}
