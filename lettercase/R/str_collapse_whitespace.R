#' Collpases multiple adjacent whitespace characters into one 
#'
#' Collapses adjacent whitespace into a single character 
#' 
#' @param string object to turn into a title case
#' @param pattern character; one or more regular expression patterns for 
#'     matching whitespace characters. Defaults to the \code{lettercase.whitespace} 
#'     option or \code{\\s, -, _} if the option has not been set.
#' @param replacement character; the whitespace character to use for replacing. 
#'     The default is "\\1" -- the character that matched.
#' 
#' @details 
#'   \code{collapse_whitespace} replaces repeated whitespace characters with a 
#'   whitespace replacement. By default, it will find and replace any of 
#'   multiple adjacent \code{pattern} characters with the replacement.  
#'   
#'   \code{pattern} can be a named patterns (see \code{?patterns}) or 
#'   character strings that are treated as regular expressions by default.
#'   
#'   To collapse mixed whitespace, provide a single patterns in the form of a 
#'   regular expression class, e.g. \code{"[\\s-_]"}. This can lead to confusing 
#'   results if a back-reference is use as a replacement value. See note below. 
#'     
#'  @note 
#'    It is inadvisable to have \code{pattern} be a multiple-character 
#'    regex class, e.g. \code{"[\\s_]"} while also using a back reference 
#'    \code{replacement} value, e.g. \code{"\\1"}. This leads to confusing 
#'    results.   
#'    
#'    \code{
#'      > str_collapse_whitespace( "A _B_ C", '[\\s-_]', "\\1" )
#'      [1] "A_B C"
#'    }                              
#'            
#'    The solution is to either provide a character vector for a \code{{pattern}}
#'    or not use back-references for the replacement depending on the desired
#'    outcome.
#'    
#'    \code{
#'      >   str_collapse_whitespace( "A _B_ C", c("\\s", "_") )
#'      [1] "A _B_ C"
#'      >   str_collapse_whitespace( "A _B_ C", '[\\s-_]', " " )
#'      [1] "A B C"
#'    }
#'    
#'                                                                                               
#' @seealso 
#'   \code{?patterns} \cr
#'   \code{\link[base]{gsub}} which is used to implement this function. \cr                                                                                                                                                                                
#'                                                                                                                                                                                                                                                                                                                                                    
#' @examples 
#'   str_collapse_whitespace( "A  B" )  
#'   str_collapse_whitespace( "A  B  C" )     
#'   str_collapse_whitespace( "A__B__C" ) 
#'   str_collapse_whitespace( "A  B__C" )   
#'   str_collapse_whitespace( "A _B_ C" )  # No transformation, no matches  
#'   
#'   # See note above:
#'   str_collapse_whitespace( "A _B_ C", '[\\s-_]' ) # possibly ill-defined 
#'   str_collapse_whitespace( "A _B_ C", c("\\s", "_") ) 
#'   str_collapse_whitespace( "A _B_ C", '[\\s-_]', " " ) 
#'
#' @export 


#? Should whitespace characters be identified in a character vector rather than
#  a string. c('\\s','-','_')  

str_collapse_whitespace <- function( 
    string
  , pattern     = getOption('lettercase.whitespace', c( pattern_whitespace, pattern_whitespace_like )  )   
  , replacement = getOption('lettercase.whitespace.replacement', "\\1" ) 
) { 
  
  for( pat in pattern ) { 
    
    pat <- paste0( '(', pat, ')' )
  
    pat <- paste0( pat, "+" )  # match multiple 
  
    string <- gsub( pat, replacement, string, perl=TRUE )
  } 
  return(string)    
}


#' @rdname str_collapse_whitespace
#' @aliases str_collapse_ws
#' @export 

str_collapse_ws <- str_collapse_whitespace

