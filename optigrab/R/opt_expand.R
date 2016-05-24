#' Expand argument vector to split names from values
#' This is an internal function and should generally not be called 
#' directly. 
#' @param opts character vector of arguments. (Default: commandArgs())
#' 
#' \code{opt_expand} does two things: 
#' \itemize{
#'   \item{"1."}{ removes values preceding and including --args }
#'   \item{"2."}{ splits and value containing an equal (=) sign }
#' }
#' 
#' @seealso 
#'   \code{\link{opt_grab}} \cr
#'   \code{\link{base}{commandArgs}}
#' 
#' @examples
#'   optigrab:::opt_expand()
#'   optigrab:::opt_expand( opts=c( "--foo", "bar") ) 
#'   
#' @note non-exported 

opt_expand <- function( opts=commandArgs() ) {
  opts <- opt_get_args(opts=opts)
  opts <- opt_split_args(opts=opts)
  
  # OPTION BUNDLING   
  # Expand dash-single-letter arguments.
  #  - splice in everything after -x as the value.
  #  wh.sd <- grep( "^-[^-]", opts ) 
  #  for( i in rev(wh.sd) ) {
  #    p <- sub( "^-[A-z]", "", opts[[i]] ) 
  #    append( opts, p, i ) 
  #  }
 
  return(opts) 
  
}
