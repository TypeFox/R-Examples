#' Getters for ANFIS object
#'
#' Obtain ANFIS's slot information, according to the given function call.
#' 
#' @param object ANFIS class object
#'
#' @return according to the call one of the following objects can be returned
#'  \item{matrix}{numeric matrix with rules or consequents}
#'  \item{list}{list with MembershipFunctions or premises and consequents}
#'  \item{character}{name of the trainingType}
#'  \item{numeric}{numeric vector with trainnig errors, fitted training values 
#'    and residuals}
#'
#' @include Anfis-initialize.R
#' @exportMethod getRules
#' @docType methods
#' @name getRules
#' @rdname ANFIS-getters
#' @aliases getRules-methods
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setGeneric(name="getRules", def=function(object){
  standardGeneric("getRules")
})
#'
#' @name getRules
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getRules,ANFIS-method
setMethod(f="getRules", signature="ANFIS", definition=function(object){
  return(object@rules)
})
#'
#' @exportMethod getPremises
#' @docType methods
#' @name getPremises
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getPremises-methods
setGeneric(name="getPremises", def=function(object){
  standardGeneric("getPremises")
})
#'
#' @name getPremises
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getPremises,ANFIS-method
setMethod(f="getPremises", signature="ANFIS", definition=function(object){
  return(object@premises)
})
#'
#' @exportMethod getConsequents
#' @docType methods
#' @name getConsequents
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getConsequents,ANFIS-method
setGeneric(name="getConsequents", def=function(object){
  standardGeneric("getConsequents")
})
#'
#' @name getConsequents
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getConsequents,ANFIS-method
setMethod(f="getConsequents", signature="ANFIS", definition=function(object){
  return(object@consequents)
})
#'
#' @exportMethod getErrors
#' @docType methods
#' @name getErrors
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getErrors,ANFIS-method
setGeneric(name="getErrors", def=function(object){
  standardGeneric("getErrors")
})
#'
#' @name getErrors
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getErrors,ANFIS-method
setMethod(f="getErrors", signature="ANFIS", definition=function(object){
  return(object@errors)
})
#'
#' @exportMethod getTrainingType
#' @docType methods
#' @name getTrainingType
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getTrainingType,ANFIS-method
setGeneric(name="getTrainingType", def=function(object){
  standardGeneric("getTrainingType")
})
#'
#' @name getTrainingType
#' @rdname ANFIS-getters
#' @inheritParams getRules
#' @aliases getTrainingType,ANFIS-method
setMethod(f="getTrainingType", signature="ANFIS", definition=function(object){
  return(object@trainingType)
})
