#' Converts a formula to character
#' 
#' Convers a formula to character representaion
#' 
#' Coerces formula to a character by deparsing.
#' 
#' @aliases as.character as.character.formula
#' 
#' @param x formula object
#' @param ... further arguments passed to or from other methods.
#'
#' @return A character vector
#' 
#' @author Christopher Brown
#' @seealso \code{\link[base]{deparse}}
#' @keywords manip utilities
#' @examples
#' 
#'   as.character( y ~ mx +  b )
#' 
#' ## The function is currently defined as
#' function(x)
#'   Reduce( paste, deparse(x) )
#' @export

as.character.formula <- function(x, ...) { 
  form <- paste( deparse(x), collapse=" " )
  form <- gsub( "\\s+", " ", form, perl=FALSE ) # remove multiple spaces
  return(form)
}
