#' str_transform
#' 
#' Convert a string by applying one or more str_* functions.
#' 
#' @param string vector. This must be an atomic vector, and will be coerced to 
#' a character vector 
#' 
#' @param ... one or more functions to apply to the string
#' 
#' @details
#'   
#'  \code{str_transform} applies successive functions to its first argument, 
#'  \code{string}.  
#'   
#' @return a character vector
#' 
#' @seealso 
#'   \code{\link{make_str_replace}}
#'   \code{\link{make_str_delete}}
#' 
#' @examples 
#'   string = c( "catch-22", "finnegans wake" )
#'   str_transform( string, str_capitalize )
#'   
#'   str_transform( string, str_capitalize, str_delete_nonword )
#'   str_delete_nonword( str_capitalize( string ) )      # SAME
#'   
#'   \dontrun{
#'     # magrittr:
#'     string %>% str_capitalize %>% str_delete_nonword   # SAME
#'   }
#'   
#' @export
           
str_transform <- function( string, ... ) { 
  
  if (!is.atomic(string)) 
    stop("String must be an atomic vector", call. = FALSE)
  if (!is.character(string)) 
    string <- as.character(string)
   
  for( fn in list(...) ) { 
    if( ! is.function(fn) ) stop( 'string transformation is not a function' )
    string <- fn( string )
  }
  
  return(string)
}
