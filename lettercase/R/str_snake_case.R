#' str_snake_case
#'
#' Function used to convert character vectors to snake case format. 
#' 
#' @param string object to turn into a title case
#' @param whitespace regular expression pattern to match for white-space
#' 
# @param acronyms character; tokens to capitalize
#' 
#' * characters are all lower case 
#' * non- \code{\\w}, \code{\\s} and - are dropped
#' * \code{\\w} and - are converted to underscore
#' * no support for acronyms
#' * multiple adjacent undescores are replaced by single underscore 
#' * Underscores at beginning or end of names are dropped  
#' 
#' @examples
#'   str_snake_case( "One Flew Over The Cuckoo's Nest" )
#'   str_snake_case( "Catch-22" )  # catch_22
#'   str_snake_case( "Catch.22" )
#'   str_snake_case( "Catch_22" )
#'   str_snake_case( "Catch  22" )
#'   str_snake_case( " Catch 22 " )
#'   
#' @rdname str_snake_case
#' @aliases str_snake_case
#' @export

str_snake_case <- function( 
    string
  , whitespace = getOption('lettercase.whitespace', '[^\\w\\s-\\.]' ) 
) { #
  
  if( ! is.character(string) ) stop( as.character(sys.call())[-1], ' is not character' )
  
  # for( ac in acronyms )  string <- gsub( tolower(ac), ac, string )
  
  string <- gsub( '[^\\w\\s-\\.]', '', string, perl=TRUE )
  string <- tolower(string)
  string <- gsub( pattern_separators, '_', string, perl=TRUE )
  string <- gsub( '__+', '_', string, perl=TRUE )  # Replace multiple with single _'s
  string <- gsub( '^_+', '', string, perl=TRUE )   # Drop leading _'s if exist
  string <- gsub( '_+$', '', string, perl=TRUE )   # Drop trailing _'s if exist
  
  return(string)
  
}

