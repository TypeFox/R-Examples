#' Unregister a an operator.
#' 
#' \code{removeOperator} unregisted an operator by removing it from the list of
#' operators found at \code{ .Options$operators }. All operator attributes are
#' that have been set will be lost.
#' 
#' Warns if the operator has not been registered.
#' 
#' @param x \code{character}. The name of the operator
#' 
#' @return None. Used for side-effects.
#' 
#' @author Christopher Brown
#' @seealso 
#'   \code{\link{setOperators}} for registering Operators.
#'   
#' @keywords utilities
#' @examples
#' 
#'   # Unregister ? as an operator.
#'   \dontrun{
#'     removeOperator( '?' )
#'   }
#' 
#' @export
removeOperator <- function(x) {

    x <- as.character(x) 
    ops <- .Options$operators 
    
    if( x %in% names(ops) ) {
      ops[[x]] <- NULL 
      options( operators = ops )
    } else {
      warning( x, " is not a REGISTERED operator.  See ?setOperator." ) 
    }

}
