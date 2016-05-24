#' Registers an operator for use with operator.tools package.
#' 
#' \code{setOperator} registers a user-defined operator as a given type.
#' Subsequently, this operator can be treated as a member of a class of
#' operators.
#' 
#' \code{setOperators} scans defined functions looking for any that have been
#' defined by the user using the special \code{any} syntax.  If found, these
#' are registered with \code{setOperator} and given the default type='user'.
#' 
#' 
#' \code{setOperator} registers a single operator similar to the way that
#' \code{setMethod} registers a method.  The definition for these operators are
#' defined by \code{.Options$operators}.
#' 
#' \code{setOperators} scans the environments for user-defined operators. If
#' found and not already registered, these are registered by
#' \code{setOperator}. Registered operators are much more efficient than
#' unregisted ones, so it is often advantageous to register the operators. When
#' \code{...} is supplied, these attributes are set for all unregistered
#' operators.
#' 
#' Operators are allowed to have attributes. The one required attribute is
#' \code{type}, which is just a character value that serves to classification
#' the operator.  On package load, All operators from base R are assigned a
#' core type as specified in \link[base]{Syntax}.  These are: namespace,
#' component, indexing, sequence, arithmetic, arithmetic, relational, logical,
#' tilde, assignment, help.
#' 
#' Users may use one of these types or assign a type of their own choosing. The
#' \code{type} is largely unrestricted, but cannot be one of the reserved
#' operator groupings: ALL, REG(ISTERED), UNREG(ISTERED), SPECIAL or user.
#' These have special meaning as described in \code{\link{operators}}. Users
#' are encouraaged to make their own types in lower case.
#' 
#' @aliases setOperator setOperators
#' @param name A character vector containing the \bold{names} of one or more
#' functions which will be registered.
#' @param type The type of operator. See Details.
#' @param ...  Attributes for the operator(s).
#' @return None. This function exists for assigning a operator to
#' \code{options('operators')}.
#' @author Christopher Brown
#' @seealso \code{\link{operators}}, \link{Syntax}
#' @keywords utilities
#' @examples
#' 
#'   \dontrun{
#'     setOperator( '%!in%', 'relational' )
#'     operators( type='relational' )
#'   }
#' @export
setOperator <- function( name, type="user", ... ) {

  if(
    type %in% 
    c( 'ALL', 'SPECIAL', 'REGISTERED', 'REG', 'UNREGISTERED', 'UNREG' ) 
  )
    stop( type, "  is a reserved operator type." )

  opts <- if( ! is.null(.Options$operators) ) .Options$operators else  list()

  # WARN IF NEW TYPE SET.
  # if( length(opts) > 1 ) {
  #  types <- sapply( opts, function(x) x[['type']] ) 
  #  if( ! type %in% types ) warning( "Adding a new operator type: ", type ) 
  # }

  # opts[[name]] <- list( name=name, type=type, fun=eval(as.name(name)) ) 
  opts[[name]] <- list( name=name, type=type ) 
  # opts[[name]][["type"]] <- type

  li <- list(...) 
  for( nm in names(li) ) opts[[name]][[nm]] <- li[[nm]] 

  invisible( base::options( operators = opts ) )
  
}

#' @export
#' @rdname setOperator
setOperators <- function( ... ) {

  ops <- apropos( "%.*%" ) 

  for( op in ops ) 
    if( ! op %in% names( .Options$operators ) ) 
      setOperator( name=op, ... ) 

}  
