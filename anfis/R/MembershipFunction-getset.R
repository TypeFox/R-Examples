#' Modify membership function parameters
#'
#' Get/set membership function parameters. 
#' 
#' @param x MembershipFunction class heirs
#' @param i numeric or character to access parameters vector [i]
#' @param value numeric parameter/s values
#'    
#' @return 
#'  \item{numeric}{parameter/s in the case of object[i]}
#'  \item{object}{MembershipFunction object in the case of object[i]<- value}
#'
#' @include MembershipFunction.R 
#' @exportMethod [
#' @docType methods
#' @name extract-methods
#' @rdname extract-methods
#' @aliases [,MembershipFunction-method
#' @seealso \code{\link{MembershipFunction-class}}
#' @family Membership Functions
setMethod(f="[",signature = "MembershipFunction", definition=function(x, i){
  return(x@parameters[i])
})
#'
#' @exportMethod [<-
#' @name extract-methods
#' @rdname extract-methods
#' @aliases [<-,MembershipFunction-method
setReplaceMethod(f="[",signature = "MembershipFunction", 
  definition = function(x,i,value){
  x@parameters[i] <- value
  validObject(x)
  return(x)
})
