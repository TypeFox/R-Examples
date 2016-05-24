utils::globalVariables('.')

#' Split a string bases on whitespace 
#' 
#' Split a string based on whitespace ignore single- and double quoted entries
#' 
#' @param x character; string to parse as if it is a command line 
#' 
#' This is an internal function used predominantly for testing. It might be 
#' deprecated in the near future.
#' 
#' @return 
#'   A character array that could be similar to that  provided by 
#'   \code{commandArgs}. 
#'   
#' @seealso 
#'   \code{\link[base]{commandArgs}}
#' 
#' @examples
#' 
#'   \dontrun{ 
#'     str <- 'cmd -t "Say Anything" --character \'Lloyd Dobler\''
#'     str_to_opts(str)
#'     split_ws_nonquote(str)
#'   }
#'   
#' @note not-exported, by design

str_to_opts <- function( x=character() )
  if( length(x) == 0 ) 
    character(0) else
    split_ws_nonquote(x)



#' @import stringi 
#' @importFrom magrittr %>%

split_ws_nonquote <- function(x) { 
  splits <- 
    "'[^']*'|\"[^\"]*\"|[^\\s]+" %>%
    stringi::stri_extract_all_regex( x, . ) 
  
  if( ! length(splits[[1]]) > 1 ) return(x)
  
  splits %>%
    magrittr::extract2(1) %>%
    stringi::stri_replace_all_regex( ., "^[\"']|[\"']", "" )
    
}
