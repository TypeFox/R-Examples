#' Split command arguments
#' 
#' Splits command argument vector to name, value pairs. 
#' This is an internal function and should generally not be called 
#' directly. 
#' 
#' @param opts character; vector of arguments. (Default: commandArgs())
#' 
#' \code{opt_split_args} splits and value containing an equal (=) sign
#' 
#' @seealso 
#'   \code{\link{opt_grab}} \cr
#'   \code{\link{base}{commandArgs}}
#' 
#' @examples
#'   optigrab:::opt_split_args()
#'   optigrab:::opt_split_args( opts=c( "--foo=hello", "-b=goodbye") ) 
#'   
#' @note non-exported 

opt_split_args <- function( opts=commandArgs()) {

  # EXPAND/Split  '=' 
  #  - Options defined with an '=', such as '-a=5' or '--alpha=5'
  #    are split into '-a' '5' and '--alpha' '5', respectively.
  
    wh.eq <- grep( "=", opts )    
    for( i in rev(wh.eq) ) {
      name.val <- strsplit( opts[[i]], "=" )[[1]]     
      opts[[i]] <- name.val[1] 
      opts <- append( opts, name.val[[2]], i ) 
    }
    return(opts)  
    
}
