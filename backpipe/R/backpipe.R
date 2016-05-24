#' backpipe
#' 
#' @description 
#' Creates backpiping operators 
#' 
#' @param pipe character; string representing the existing pipe operator
#' @param backpipe character; string representing the desired backpipe operator
#' 
#' @details
#' Only \code{pipe} is necessary. Arbitrary mixing of forward and reverse are 
#' not allowed.
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
#'   backpipe('%>%')
#'   backpipe('%>>%', '%<<%')
#'                     
# @rdname 
# @aliases %<% %<<% 
#' @export

backpipe <- function(pipe, backpipe=gsub(">","<",pipe) ) { 

  function( lhs, rhs ) { 
  
    lhs <- substitute(lhs)
    rhs <- substitute(rhs)
    
    parent = parent.frame() 
      
    # Allow for nested backpipes 
    #  The first case recomposes the nested backpipe expressions, 
    #  The second part uses forward pipe for evaluation to keep this close to the 
    #  original
    ca <- 
      if( is.call(lhs) && deparse( lhs[[1]] ) == backpipe ) {
        rhs.  <- call( pipe, rhs, lhs[[ length(lhs)]] )
        lhs. <- lhs[[2]]
      
        call( backpipe, lhs[[2]], rhs. )
  
      } else { 
        call( pipe, rhs, lhs )
      }
    
    eval(ca, parent, parent )
  }
}
