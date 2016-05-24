#' Get verb from the command line
#' 
#' Return the verb for the application. The verb is the first argument that is 
#' not part of an option. 
#' 
#' @param opts character; Vector from which to parse options 
#'        (default: \code{commandArgs()} )
#'
#' @details
#' 
#' Some applications such as \emph{git} support commend verbs, e.g. \code{push},
#' \code{fetch}, etc. These can be retrieved by \code{opt_get_verb}.
#' 
#' \code{opt_get_verb} look for the first unaccounted for options (after 
#' \code{--args}). The number of options needed by each flags is determined by 
#' and saved by calls to \code{\link{opt_get}}. See the details to see how 
#' \code{n} is automatically determined.  When not explicitly defined the number
#' of options needed by each flag is 1.  Becasue of this, it convention to call 
#' \code{opt_get_verb} after all \code{opt_get} calls.  For most simple 
#' applications, it likely doesn't matter.  
#' 
#' @section Assumptions:
#' 
#' \code{opt_get_verb} assumes any flags occurring before the verb have 
#' exactly 1 value. A command line such as "> myscript --verbose verb" will be 
#' misparsed; the code will assume that "verb" is the value of the 
#' flag "--verbose"
#'  
#' @return character of length 1; the verb found from the command-line. 
#' \code{NA} if a verb cannot be identified.
#' 
#' @seealso 
#'   \code{\link{opt_get}} \cr
#'   \code{\link{base}{commandArgs}}
#' 
#' @examples
#'   opt_get_verb()
#'   
#' @export

opt_get_verb <- function( opts=commandArgs() ) {
  
  # browser()
  opts <- opt_get_args(opts=opts)
  
  if ( length(opts) == 0 ) return(NA)
  
  
  found_flag <- FALSE
  
#   wh <- opts %>% is.flag %>% which %>% tail(1) 
#  
#  # No FLAGS FOUND ... first option must be a command
#  if( length(wh) == 0 ) return( opts[[1]])
#  
#  flag <- opts[[wh]]
# 
#  # How many arguments are expected from `flag`
#  ops <- getOption('optigrab')$options
#  if( exists(flag, ops) ) 
#    n <- ops[[flag]][['n']] else 
#    n <- 1
  
#  wh.command <- wh + n + 1
#  if( wh.command > length(opts) ) return(NA)
  
  # return( )
  
  ops <- getOption('optigrab')$options
#   i <- 1
#   while( i <= length(opts) ) { 
#     if( is.flag(opts) )
#     ops  
#     
#   }
  

  n = 0     # how many options needed, starts at 0.  
  for(opt in opts) {
    
    if( is.flag(opt) ) { 
      if( exists( opt, ops ) )
      n <- ops[[opt]][['n']] else 
      n <- 1
      next      
    } 
    
    if( n==0 )          
      return(opt) else  # no options needed
      n <- n -1         # give-up on option
  }  
  
  return(NA_character_) # no verb found
  
}
