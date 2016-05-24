# @param x character; vairaible name or flag to transform or test.

java_flag_test <- function(x) grepl( '^-\\S*$', x )

java_flag_to_name <- function(x) gsub( "^-", "", x)

java_name_to_flag <- function(x)  paste0( "-" , x ) 


#' Java-style command line options
#' 
#' Functions for enabling Java-style commpand-line option behavior. 
#'  
#' @details 
#'   Functions for enabling Java-style command-line options. Java-style options
#'   are characterized by a single dash (\code{-}) before the option name.  
#'      
#' @seealso 
#'   Non-exported function \code{*_flag_test}, \code{*_flag_to_name} and 
#'   \code{*_name_to_flag} \cr
#'   \code{\link{gnu_style}} \cr
#'   \code{\link{java_style}} \cr
#'   \code{\link{ms_style}}
#'      
#' @rdname java_style
#' @export

  java_style = list( 
      flag_test        = java_flag_test
    , flag_to_name     = java_flag_to_name 
    , name_to_flag     = java_name_to_flag
  )
