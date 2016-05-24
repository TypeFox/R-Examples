#' backpipes: \%<\%, \%<<\%
#' 
#' @description 
#' Provides back-piping operators \code{\%<\%} (magrittr) or 
#' \code{\%<<\%} (pipeR) allowing a reverses (right-to-left) order of 
#' arguments
#' 
#' @param lhs argument on the lhs of the backpipe
#' @param rhs argument on the rhs of the backpipe
#'
#' @details
#' 
#' * \code{\%<\%} works with magrittr. \code{\%<<\%} works with piper. They are
#'   both based on \code{\link{backpipe}}.\cr
#' 
#' * Other magrittr operators and complex expressions are not supported yet. 
#' 
#' * It is not possible to mix forward and and backward piping in the same 
#'   expression because of likely ambiguous results.
#' 
#' @references 
#'   \url{https://github.com/smbache/magrittr/issues/26} \cr
#'   \url{http://stackoverflow.com/questions/31305342/is-right-to-left-operator-associativity-in-r-possible} \cr
#'   
#' @seealso 
#'   \code{\link[magrittr]{\%>\%}} \cr
#'   \code{\link[pipeR]{\%>>\%}} \cr
#'   
#' @examples
#'   \dontrun{
#'     require(magrittr)
#'     letters %>% paste0( 1:26 )  # forward pipe
#'     paste0( 1:26 ) %<% letters  # backward pipe
#' 
#'     mean %<% range( na.rm = TRUE ) %<% c(1:5, NA)
#'   }
#'   
#'   \dontrun{
#'     require(pipeR)
#'     letters %>>% paste0( 1:26 )  # forward pipe
#'     paste0( 1:26 ) %<<% letters  # backward pipe
#' 
#'     mean %<<% range( na.rm = TRUE ) %<<% c(1:5, NA)
#'   }
#'   
#'   \dontrun{
#'     require(shiny)
#'     div( class="outer-outer") %<%
#'       div( class="outer")     %<% 
#'         div( class="inner")   %<% 
#'           h1( "content", role="heading" )
#'    }          
#'
#' @export 
#' @rdname backpipes
#' @aliases %<%                       
`%<%` <- backpipe('%>%')


#' @rdname backpipes
#' @aliases %<<%
`%<<%` <- backpipe('%>>%')



# This is conditional export based on whether magrittr or pipeR are loaded.
# Support magrittr: It also requires:
#     exportPattern( "%.*%")
# in the NAMESPACE file
#  There may be other magrittr packages
# Conditional assignment
# if( any( grepl( "^package:magrittr.*", search() ) ) ) 
#  assign( '%<%', backpipe('%>%') )

# Support pipeR
# if( "package:pipeR" %in% search() ) 
#   assign( '%<<%', backpipe('%>>%') )
