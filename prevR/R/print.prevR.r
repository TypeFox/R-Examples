#' Summary of a prevR object.
#' 
#' Method \code{print} for objects of class \code{\link[=prevR-class]{prevR}}: 
#' shows a summary of the object's characteristics.
#' 
#' @param x object of class \code{\link[=prevR-class]{prevR}}.
#' 
#' @note Exactly the same as \code{\link{show,prevR-method}}.
#' @seealso \code{\link{summary,prevR-method}}.
#' @examples 
#' print(fdhs)
#' \dontrun{
#'  dhs <- rings(fdhs,N=c(100,300,500))
#'  print(dhs)
#' }
#'
#' @aliases print print-methods print,prevR-method

setMethod("print","prevR",
function(x){
   show(x)
   invisible(NULL)
  }
)
