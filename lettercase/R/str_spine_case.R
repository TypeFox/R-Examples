#' str_spine_case
#'
#' Function used to convert character vectors to spine case format. 
#' 
#' @param string object to turn into a title case
#' @param whitespace character; regular expression pattern for matching whitespace
#'   characters
#' @param collapse.whitespace logical; whether adjacent whitespace is collapsed
#' 
# @param acronyms character; tokens to capitalize
#' 
#' * characters are all lower case 
#' * non- \code{\\w}, \code{\\s} and - are dropped
#' * \code{\\w} and - are converted to underscore
#' * no support for acronyms
#' 
#' @examples
#'   str_spine_case( "One Flew Over The Cuckoo's Nest" ) # One Flew Over The Cuckoo's Nest
#'   str_spine_case( "Catch 22" )  # catch-22
#'   str_spine_case( "Catch_ 22" )  
#'   
#' @rdname str_spine_case
#' @aliases str_spine_case
#' @export

  str_spine_case <- 
    function( string
            , whitespace          = getOption('lettercase.whitespace', '[\\s-_]') 
            , collapse.whitespace = getOption('collapse.whitespace', TRUE )
    ) { 
    
      if( ! is.character(string) ) stop( as.character(sys.call())[-1], ' is not character' )
      
      # for( ac in acronyms )  string <- gsub( tolower(ac), ac, string )
      
      string <- gsub( '[^\\w\\s-]', '', string, perl=TRUE )
      
      if( collapse.whitespace ) whitespace <- paste0( whitespace, "+" )
      string <- gsub( whitespace, '-', string, perl=TRUE )
      
      return(string)
    
  }


#' @rdname str_spine_case
#' @aliases str_spinal_case
#' @export 

  str_spinal_case <- str_spine_case


#' @rdname str_spine_case
#' @aliases str_hyphen_case
#' @export 

  str_hyphen_case <- str_spine_case
