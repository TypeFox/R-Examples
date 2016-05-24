#' \code{Print} a MembershipFunction object
#'
#' Generic Print Method for MembershipFunction class and descendants. 
#' Usage: print(x, ...)
#'
#' @param x MembershipFunction class object
#' @param ... not used but included for generic print compatibility
#'
#' @return console output of the object
#'
#' @include MembershipFunction-show.R
#' @exportMethod print
#' @docType methods
#' @aliases print,MembershipFunction-method
#' @seealso \code{\link{MembershipFunction-class}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setMethod(f="print",signature="MembershipFunction",
  definition = function(x,...){
    show(x)
  }
)