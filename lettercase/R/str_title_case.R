#' str_title_case
#'
#' internal function for converting character to title case.  This is not 
#' exported and should not be called directly. 
#' 
#' @param x object to turn into a title case
# @param acronyms character; tokens to capitalize
#' 
#' @examples
#'   str_title_case( "One Flew Over The Cuckoo's Nest" )
#'   str_title_case( "one_flew_over_the_cuckoo_'_s_nest" )
#'   
#' @aliases str_title_case
#' @export

str_title_case <- function(x ) {
  
  # for( ac in acronyms )  x <- gsub( tolower(ac), ac, x )
  
  x <- gsub( "[\\s_]+", " ", x, perl=TRUE )  # whitespace to single space
  
  x <- gsub( "\\b([a-z])([a-z]+)", "\\U\\1\\L\\2" ,x, perl=TRUE ) # ucfirst
  
  return(x)
  
}
