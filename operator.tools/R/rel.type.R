#' Get the relational type of a relational operator.
#' 
#' \code{rel.type} gets the relational type of a relational operator. The
#' relational type is one of \code{'gt'}, \code{'lt'}, \code{'eq'},
#' \code{'ne'}.
#' 
#' A relational operator is an operate that relates the relationship between
#' arguments. The core relational operators are: >, >=, <, <=, ==, !=, %in%.
#' 
#' The relational.type is a simple roll-up of these operators. > and >= are gt,
#' etc. The value is retrieved from .Options$operators[[x]][['rel.type']] and
#' can be defined for relational operators using \code{\link{setOperator}}.
#' 
#' A relational type provides an indication of nature of the relational
#' operator.
#' 
#' @aliases 
#'   rel.type rel.type.function rel.type.name rel.type.call
#'   rel.type.expression
#' 
#' @param x An operators expressed as a \code{function} or \code{name}
#' 
#' @return \code{character} value of the operator.  One of: \code{'gt'},
#' \code{'lt'}, \code{'eq'}, \code{'ne'}.
#' 
#' @author Christopher Brown
#' 
#' @seealso 
#'   \code{\link{operators}}, \code{\link{setOperator}}
#'   
#' @keywords utilities
#' 
#' @examples
#'  \dontrun{
#'   rel.type( `==` ) 
#'   rel.type( as.name('==') )
#'  }
#' 

#' @export 
rel.type <- function(x) 
  UseMethod( 'rel.type', x ) 

#' @export 
rel.type.name <- function(x) 
  if( is.operator(x, type='relational') ) {
    .Options$operators[[as.character(x)]][['rel.type']]
  } else { 
    NULL
  }

#' @export 
rel.type.function <- function(x) 
  if( is.operator(x, type='relational') ) {
    .Options$operators[[ fun2name(x) ]][['rel.type']] 
  } else { 
    NULL
  }

#' @export 
rel.type.call <- function(x) 
  if( is.name(x[[1]]) ) 
    rel.type( as.name( x[[1]] )  ) else
    NULL

#' @export   
rel.type.expression <- function(x) 
  sapply(x,rel.type)  


# rel.type expression.  This does not quite work since expressions
# are not necessarily a LHS OP RHS construction. 
# rel.type.expression <- function(x) 
#   sapply(x, function(x) rel.type(x[[1]])  )

