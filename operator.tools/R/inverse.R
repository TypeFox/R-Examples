#' Invert an R operator
#' 
#' \code{inverse} is a S3 generic method for inverting an R operator in the
#' mathematical sense. Presently, inverses are defined for relational
#' operators, i.e. changing \code{>} to \code{<=} etc.
#' 
#' Arguments will be checked against the defined list of inverses,
#' If an entry exists, the corresponding inverse is returned.
#' 
#' @name inverse
#' 
#' @aliases inverse inverse inverse.name inverse.function
#' 
#' @param x object representing an R operator
#' @param ... additional arguments 
#' 
#' @return 
#'   \code{inverse} returns the inverse in the same form as the \code{x}
#'   argument. Thus, if a name is provided, a name is returned. If a function is
#'   provided, a function is returned.
#'   
#' @author Christopher Brown
#' 
#' @seealso 
#'   \code{\link{operators}} especially \code{operators(type="relational"))}
#'   
#' @references http://en.wikipedia.org/wiki/Inverse_mathematics.
#' 
#' @keywords methods symbolmath utilities
#' 
#' @examples
#'   \dontrun{
#'     inverse( as.name( '!=' ) )
#'     inverse( `==` )
#'  } 

#' @export
inverse <- function(x, ... ) 
  UseMethod( 'inverse' ) 

#' @export
inverse.name <- function(x, ...) {
  
  op <- as.character(x)
  inverses <- sapply( .Options$operators, function(x) x$inverse )
  # inverses <- inverses[ ! sapply( inverses, is.null ) ] 

  ret <- inverses[[ op ]]
  if( ! is.null(ret) ) return( as.name(ret) ) 

  # if( op %in% names(inverses) ) {
 #   return( as.name( inverses[[ op ]] ))
 # } else  {
  warning( "No inverse found for op: ", op, call.=FALSE )
  # }

  return(ret)
}


#' @export
inverse.function <- function(x, ...) {

  inverses <- sapply( .Options$operators, function(x) x$inverse )

  return( name2fun( as.name(inverses[[fun2name(x)]]) ) )
  warning( "No operator matched" )
  return( NULL )

}
