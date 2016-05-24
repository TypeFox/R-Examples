#' str_delete
#' 
#' Delete or remove characters from a string based on one or more patterns
#' 
#' @param string atomic character vector 
#' @param ... stringr-style matching \code{\link[stringr]{modifiers}} 
#' 
#' @details 
#' 
#' Deletes all occurences of the patterns from the string using 
#' \code{\link[stringr]{str_replace_all}}
#' 
#' 
#' @references 
#'  \code{\link[stringr]{modifiers}} \cr
#'  \code{\link[stringr]{str_replace_all}}
#' 
#' @examples
#'   
#'   str_delete( "ABC & 123", stringr::regex("\\W") )  # ABC123
#' 
#' @import stringr
#' @export 

str_delete <- function( string, ... ) { 

  for( pattern in list(...) ) {
    string <- stringr::str_replace_all(string,pattern,'')
  }
  return(string)    

}
