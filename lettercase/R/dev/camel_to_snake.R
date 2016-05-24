#' camel_to_snake
#' 
#' Change CamelCase to snake_case
#'
#' @param x character; vector to change from CamelCase 
#' @param acronyms logical; whether to detect and preserve acronyms case 
#' (default: \code{FALSE})
#' 
#' Changes CamelCase character vectors to snake_case.  The trna e.g FooBar  Treats consecutive 
#' capital letters as an aabreviation and does not coer
#' 
#' @author Christopher Brown
#' 
#' @note copyright (c)2013 used with permission 
#' 
#' @references
#'   http://titlecase.com/
#'   http://en.wikipedia.com/wiki/Camel_case
#'   http://en.wikipedia.com/wiki/Snake_case
#'   
#' @seealso 
#'   \code{ \link{make_names} } 
#'     
#' @examples
#'   
#'   camel_to_snake( "FooBar" )     # foo_bar
#'   camel_to_snake( "fooBar" )     # foo_bar
#'   camel_to_snake( "fooBAR" )     # foo_BAR  
#'   camel_to_snake( "FooBarBaz" )  # foo_bar_baz
#'   camel_to_snake( "fooBarBaz" )  # foo_bar_baz 
#'   
#'   # Embedded Acronyms 
#'   camel_to_snake( "fooBARBaz", acronyms=TRUE )  # foo_BAR_baz 
#'   camel_to_snake( "fooBARBaz", acronyms=TRUE  )  # foo_BAR_baz
#'   camel_to_snake( "FOOBarBaz", acronyms=TRUE  )  # foo_bar_baz
#'    
#'   camel_to_snake( "foo_bar_baz" ) # foo_bar_baz 
#' @export

camel_to_snake <- function( x, acronyms = FALSE ) { 
  
  # 4 gsubs
  #  1. First position: Ax -> ax 
  #  2. aBx -> a_bx
  #  3. aBC -> a_BC 
  #  4. ABc -> A_bc
  
  x <- gsub( '^([A-Z])([^A-Z])'      , '\\L\\1\\2'    , x, perl=TRUE )
  x <- gsub( '([a-z])([A-Z])([^A-Z])', '\\L\\1_\\2\\3', x, perl=TRUE )
  x <- gsub( '([a-z])([A-Z])([A-Z])' , '\\1_\\2\\3'   , x, perl=TRUE )
  x <- gsub( '([A-Z])([A-Z])([a-z])' , '\\1_\\L\\2\\3', x, perl=TRUE )
  
  if( ! acronyms ) x <- tolower(x)
  
  return(x)
  
}  
