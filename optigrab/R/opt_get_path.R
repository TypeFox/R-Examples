#' Get path current running script
#' 
#' @param opts character; cector from which to parse options 
#'   (default: \code{commandArgs()} )
#'   
#' @param full.name boolean; expand to full path(?)
#'  
#' \strong{ This function is deprecated, use \code{this_file} instead.}   
#'     
#' @return character; path to Rscript or \code{NA} if there isn't one. 
#' 
#' @references
#'   \url{http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script}
#' 
#' @seealso 
#'   \code{\link{opt_grab}} \cr
#'   \code{\link{base}{commandArgs}}
#' 
#' @examples
#'   optigrab:::opt_get_path()
#'   
#' @export

opt_get_path <- function( opts=commandArgs(), full.name = FALSE ) {
  
  warning( "'opt_get_path' is deprecated. Use 'this_file' instead.")
  return( this_file( opts, full.path = full.name) )
  opts <- opt_split_args(opts)
  
  wh.args <- grep( "--file", opts )[1]  # i.e. first occurence of --file
  
  if ( is.na(wh.args) || (wh.args == length(opts)) ) return(NA)

  path <- opts[ wh.args + 1 ]
  if( full.name ) 
    return( normalizePath( path, mustWork = FALSE ) ) else
    return( path )
  
}
