#' str_ucfirst
#' 
#' Convert first character of a string to uppercase
#' 
#' @param string character vector to be converted.
#' 
#' @seealso 
#'   \code{\link{str_title_case}}
#'   
#' @examples
#'  str_ucfirst( "one flew over the cuckoo's nest" )
#'  str_ucfirst( "catch-22" )  
#'  str_ucfirst( "Portrait of the Artist as a Young Man" )
#'  
#' @export

str_ucfirst <- function(string) {
  
  if( ! is.character(string) ) 
    stop( deparse(substitute(string)), ' is not character' )
  
  gsub( "\\b([a-z])([a-z]+)", "\\U\\1\\L\\2", string, perl=TRUE ) # ucfirst

}


# @rdname string_transformations
# @export
#    str_ucfirst <- make_str_replace( "\\b([a-z])([a-z]+)", "\\U\\1\\L\\2" )


#' @rdname str_ucfirst
#' @export
#' @aliases is_ucfirst 
#' @examples 
#' # is_ucfirst 
#'   is_ucfirst( 'ABC123' )      # TRUE
#'   is_ucfirst( 'abc123' )      # FALSE 
#'   is_ucfirst( 'aBC'  )        # FALSE
#'   is_ucfirst( 'Abc' )         # TRUE 
#'   is_ucfirst( 'Abc dEF' )     # FALSE 
#'   is_ucfirst( '123' )         # TRUE 
#'     
is_ucfirst <- function(string) 
  ! grepl( "\\b([a-z])", string )
