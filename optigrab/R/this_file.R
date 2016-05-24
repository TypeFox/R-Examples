#' this_file 
#' 
#' Name or path to the current file 
#' 
#' @param opts character; vector of arguments. (Default: \code{commandArgs()})
#' 
#' @param local logical; if \code{TRUE} returns the most currently sourced 
#'    script as opposed to the orignal/first source script. 
#'    (Default: \code{TRUE}) 
#'    
#' @param full.path logical; Whetther to return the full path to the sourced 
#'    file. (Default: \code{TRUE})
#'
#' @details 
#' \code{this_file} returns the name or path of the executing file whehter 
#' the file was invoked from \strong{Rscript} or in an interactive session. 
#' Further it \code{source} 
#' 
#' Argument \code{local} controls whether it is the current file (\code{TRUE}) 
#' or the orignal, top-level file.
#' 
#' 
#' @references
#'   \url{http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script}
#'
#' @return one-element character vector with the path to the current file; 
#' returns \code{NA} is in an interactive session not in a file.
#' 
#' @seealso 
#'   \code{\link{opt_grab}} \cr
#'   \code{\link{base}{commandArgs}}
#' 
#' @examples
#'   this_file()
#' 
#' @export

this_file <- function( 
    opts      = commandArgs()
  , local     = TRUE
  , full.path = TRUE 
) { 
  # browser()
  current_script <- 
    if(local)
      sys.frame(1)$ofile
      # tail( unlist(lapply(sys.frames(), function(env) env$ofile)), 1 ) else
      NULL
  
  if( is.null(current_script) ) {   # No source, RScript?

    opts <- opt_split_args(opts)   
    wh.args <- grep( "--file", opts )[1]  # i.e. first occurence of --filele
    
    if ( is.na(wh.args) || (wh.args == length(opts)) ) return(NA)
    path <- opts[ wh.args + 1 ]
    
    if( full.path ) 
      return( normalizePath( path, winslash=.Platform$file.sep, mustWork = TRUE  ) ) else
      return( path )
      match <- grep( "--file", opts )
   
  } # else {
  
  # 'source'd via R console
  if(full.path)
    return( normalizePath( current_script, winslash=.Platform$file.sep, mustWork=TRUE ) ) else
    return(current_script)  
  
}