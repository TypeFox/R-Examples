#' Return the type for an operator.
#' 
#' Given an operator or its name/symbol, return the \bold{type} of operator.
#' 
#' The operator is first checked against all operators that have been
#' registered with the \code{\link{setOperator}} command.  If there is a match,
#' its type is returned. If no matching operator is found, \code{op} is matched
#' against unregistered operators that have been defined with the
#' \code{\%any\%}-syntax. If a match is found, UNREGISTERED is returned.
#' 
#' The list of operators are maintained in \code{.Options\$operators} and be
#' altered suing the \code{\link{setOperator}} command.
#' 
#' @aliases operator.type operator.type.name operator.type.function
#' @param op An operator either as a name/symbol or function.
#' @return A \code{character} value.
#' 
#' For registered operators, the registered type is returned.  For Base R
#' operators, the types come from \code{\link[base]{Syntax}}.
#' 
#' For operators defined with the \code{\%any\%}-syntax but, not registered using
#' \code{\link{setOperator}}, "UNREGISTERED" is returned.
#' 
#' NULL is returned otherwise.
#' @author Christopher Brown
#' @seealso \code{\link{operators}}, \code{\link{setOperator}}. \link{Syntax}
#' @keywords utilities manip methods
#' @examples
#' 
#' 
#'  \dontrun{
#'   operator.type( `+` )
#'   operator.type( `<=` )
#'   
#'   e <- quote( A +B )
#'   operator.type( e[[1]] )
#' 
#'   operator.type( as.name('+') )
#'  }
#'  


#' @export
operator.type <- function(op) 
  UseMethod( 'operator.type', op )

#' @export
operator.type.name <- function(op) {
  op <- as.character( deparse( op ) ) 

  if( op %in% operators(types="REG") )
    return( .Options$operators[[op]]$type )
 
  if( op %in% operators(types="UNREG") )
    return( "UNREGISTERED" ) 

  return(NULL)
}

#' @export
operator.type.function <- function(op) {

  # REGISTERED OPERATORS
    li.fun <- sapply( operators( types="REG"), function(x) eval(as.name(x) ) )

    for( nm in names(li.fun) ) 
      if( identical( op, li.fun[[nm]] ) ) return( operator.type( as.name(nm)) )


  # UNREGISTERED OPERATORS
    li.fun <- sapply( operators( types="UNREG"), function(x) eval(as.name(x) ) )

    for( nm in names(li.fun) ) 
      if( identical( op, li.fun[[nm]] ) ) return( "UNREGISTERED" )

    return(NULL) 

}


