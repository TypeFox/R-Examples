#' \code{Show} a MembershipFunction object
#'
#' Generic display method for MembershipFunction class and its descendants. 
#' Usage: show(object)
#' 
#' @param object MembershipFunction class object
#'
#' @return console output of the object
#' 
#' @include MembershipFunction.R
#' @exportMethod show
#' @docType methods
#' @aliases show,MembershipFunction-method
#' @seealso \code{\link{MembershipFunction-class}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setMethod(f="show",signature="MembershipFunction",
  definition = function(object){
    if(class(object)=="MembershipFunction"){
      stop("Can't display Virtual Class")
    }else{
      if(length(object@parameters)==0 | object@name==""){
        stop("No parameters or name available")
      }else{
        cat("MembershipFunction: ", object@name, "\n")
        cat("Number of parameters:", object@nParameters,"\n")
        print(object@parameters)
        cat("Expression: ");print(object@expression)
      }
    }
  }
)