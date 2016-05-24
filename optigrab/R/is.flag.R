#' Determine if/which vector element are options flags
#' 
#' Determines if an element of a vector is an option flag (as opposed to a 
#' value ) by checking against the option style
#' 
#' @details \code{is.flag} and \code{which.flag} are internal functions not 
#' expected to be called directly. 
#'  
#' They are  used to identify which elements of the option vector are 
#' option names (as opposed to option values). Options are identified by  
#' \code{getOptions('optigrab')$style$flag_test}. By defailt, \emph{optigrab} 
#' follows GNU-style command line arguments, i.e. those beginning with 
#' "--" or "-" and are set at the time of package loading.
#' 
#' @param x vector of options, for example \code{commandArgs()}.
#' @return logical. indicating which arguments are flags.
#' @examples
#'   optigrab:::is.flag( c( "--foo", "bar") )
#'   optigrab:::is.flag( c( "--foo", "bar", "-f", "-b", "text" ) )
#'   
 
is.flag <- function(x) {
  return( getOption("optigrab")$style$flag_test(x) )
}

#' @rdname is.flag
#' @return numeric
#' @examples
#'   optigrab:::which.flag( c( "--foo", "bar") )
#'   optigrab:::which.flag( c( "--foo", "bar", "-f", "-b", "text" ))
#' 

which.flag <-function(x) which( is.flag(x) )
