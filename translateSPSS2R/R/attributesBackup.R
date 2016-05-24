#' Stores Attributes 
#'
#' Attribut backup function
#'
#' The conditions to select cases are specified in a logical expression. These logical expressions can contain relational operators, logical operators and arithmetic operations. 
#' 
#'  \strong{NOTE:} For temporary case selection, specify a TEMPORARY command before SELECT IF.
#'
#' @param x a (non-empty) data.frame, data.table object or input data of class \code{"xpssFrame"}. 
#' @return Output is a subset of the actual dataset under the condiftion of the logical expression.
#' @seealso \code{\link{attributes}} \code{\link{attr}}
#' @author Bastian Wiessner
#' #' @examples
#' #load data
#' data(fromXPSS)
#' attributes(fromXPSS)
#' attributes(fromXPSS$V7_2)
#' x <- attributesBackup(fromXPSS)
#' fromXPSS <- fromXPSS[order(fromXPSS$V2),] 
#' attributes(fromXPSS)
#' attributes(fromXPSS$V7_2)
#' fromXPSS <- applyAttributes(fromXPSS, x)
#' attributes(fromXPSS)
#' attributes(fromXPSS$V7_2)
#' @export
attributesBackup <- function(x) {
  
     
  attrs <- list("global","local")
  
  attrs$global <- attributes(x)
  
  
  attrs$local <- sapply(x, function(x){
      attributes(x)
  })
  
    
  class(attrs) <- c("xpssAttributes", "list")
    
  return(attrs)
}




