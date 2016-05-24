#' str_cap_words
#'
#' Function used to convert character vectors to CapWord format. 
#' 
#' @param string character vector  to turn into a CapWords
#' 
# @param acronyms character; tokens to capitalize
#' 
#' CapWords is distinguished by:
#' \enumerate{  
#'   \item All words with upper case first character 
#'   \item No non-word characters allowed (letters/numbers only allowed)
#' }
#' 
#' The recipe for chaging into capwords is:
#' \enumerate{
#'   \item replace whitespace and other word separators with space
#'   \item str_ucfirst
#'   \item strstdelete whitespace
#' }     
#'   
#'    
#' @examples
#'   # CAP WORDS 
#'     str_cap_words( "One Flew Over The Cuckoo's Nest" )
#'     str_cap_words( "Catch-22" )  # CATCH
#'  
#' @rdname str_cap_words
#' @aliases cap_words
#' @export

str_cap_words <- function(string) { #
  
  if( ! is.character(string) ) stop( deparse(substitute(x)), ' is not character' )
  
  results <- str_delete_space( str_ucfirst( str_replace( string, pattern_separators, ' ' ) ) )
  
# results <- string %>% 
#     str_replace( pattern_separators, ' ' ) %>% 
#     str_ucfirst %>% 
#     str_delete_space
  
  return(results)
  
}


#' @rdname str_cap_words
#' @examples
#'    is_cap_words( "AbcDef" )  # TRUE 
#'    is_cap_words( "Abc" )     # TRUE 
#'    
#'    is_cap_words( "abcdef" )  # FALSE 
#'    is_cap_words( "123Abc" )  # FALSE
#'    is_cap_words( "_Abc" )    # FALSE
#'    is_cap_words( "Abc_Def" ) # FALSE 
#'    is_cap_words( "Abc$Def" ) # FALSE 
#'    
#' @export 

# is_cap_words <- function( string ) 
#   ! ( 
#     grepl( '^[^A-Z]', string, perl=TRUE ) |    # initial non-cap
#     grepl( '[_\\W]' , string, perl=TRUE )      # non-word or underscore 
#   )