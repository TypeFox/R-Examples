#' MembershipFunction S4 class
#'
#' Represent a generic virtual S4 MembershipFunction class, for fuzzy further 
#' redefinition. The actual subclases available are GaussianMF, 
#' NormalizedGaussianMF and BellMF.
#'
#' @section Functions:
#' MembershipFunction S4 class includes the following functions:
#' \describe{
#'  \item{show/print}{generic output of the object.}
#'  \item{"[", "[<-"}{getter and setter of the parameters values.}
#'  \item{evaluateMF}{return membership value at x.}
#'  \item{derivateMF}{return the derivate membership at x.}
#' }
#' 
#' @note validity: nParameters == length(parameters) and parameters != NA and 
#' names(parameters)!="".
#'
#' @slot parameters named numeric vector with parameters of Membership Function.
#' @slot nParameters integer with the number of parameters for validity check.
#' @slot name character The description of the membership function.
#' @slot expression expression object just to display purposes.
#' 
#' @name MembershipFunction-class
#' @rdname MembershipFunction-class
#' @import methods
#' @exportClass MembershipFunction
#' @seealso \code{\link{BellMF-class}}, \code{\link{GaussianMF-class}} or 
#'  \code{\link{NormalizedGaussianMF-class}}
#' @family Membership Functions
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
MembershipFunction<-setClass(Class="MembershipFunction",
  slots=list(parameters="numeric",
    nParameters="integer",
    name="character",
    expression="expression"),
  contains="VIRTUAL",
  prototype=prototype(parameters=numeric(length=0),nParameters=0L, name=""),
  validity=function(object){
    if(object@nParameters != length(object@parameters)){
      stop("nParameters != length(parameters)")
    }
    if(any(is.na(object@parameters))){
      stop("Assigning NA to parameters")
    }
    if(any(names(object@parameters)=="")){
      stop("Complete the names of the parameters")
    }
    return(TRUE)
  }
)
