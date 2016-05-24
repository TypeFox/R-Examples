#' ANFIS training results
#'
#' Obtain ANFIS slot information, according to training output
#' 
#' @param object ANFIS class object
#' @param ... required by resid, residuals, coef and coefficients
#'
#' @return according to the call one of the following objects can be returned
#'  \item{list}{list with premises and consequents.}
#'  \item{numeric}{numeric vector with training errors, fitted training values 
#'    and residuals.}
#'  \item{printed}{statistics of the training process.}
#'
#' @include Anfis-printshow.R
#' @exportMethod fitted
#' @docType methods
#' @name fitted
#' @rdname ANFIS-metrics
#' @inheritParams getRules
#' @aliases fitted,ANFIS-method
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setMethod(f="fitted", signature="ANFIS", definition=function(object, ...){
  return(object@fitted.values) 
})
#'
#' @exportMethod fitted.values
#' @name fitted.values
#' @rdname ANFIS-metrics
#' @inheritParams fitted
#' @aliases fitted.values,ANFIS-method
setMethod(f="fitted.values", signature="ANFIS", 
  definition=function(object, ...){
    return(fitted(object))
})
#'
#' @exportMethod coef
#' @name coef
#' @rdname ANFIS-metrics
#' @inheritParams fitted
#' @aliases coef,ANFIS-method
setMethod(f="coef", signature="ANFIS", definition=function(object, ...){
  return(list(premises=object@premises, consequents=object@consequents))
})
#'
#' @exportMethod coefficients
#' @name coefficients
#' @rdname ANFIS-metrics
#' @inheritParams fitted
#' @aliases coefficients,ANFIS-method
setMethod(f="coefficients", signature="ANFIS", 
  definition=function(object, ...){
    return(coef(object))
})
#'
#' @exportMethod resid
#' @name resid
#' @rdname ANFIS-metrics
#' @inheritParams fitted
#' @aliases resid,ANFIS-method
setMethod(f="resid", signature="ANFIS", definition=function(object, ...){
  return(object@residuals)
})
#'
#' @exportMethod residuals
#' @name residuals
#' @rdname ANFIS-metrics
#' @inheritParams fitted
#' @aliases residuals,ANFIS-method
setMethod(f="residuals", signature="ANFIS", definition=function(object, ...){
  return(resid(object))
})
#'
#' @exportMethod summary
#' @name summary
#' @rdname ANFIS-metrics
#' @inheritParams fitted
#' @aliases summary,ANFIS-method
setMethod(f="summary", signature="ANFIS", definition=function(object, ...){
  ##Check training
  if(object@trainingType == "Not trained yet"){
    warning("Traying to summary a not trained ANFIS")
  }else{
    show(object)
    cat("\nCall: ")
    print(object@call)
    if(object@trainingType == "trainHybridJangOnLine"){
      cat("\nStatistics for On-line training\n")
      cat("\nPattern errors\n")
      print(summary(object@errors))
      cat("\nEpoch errors\n")
      epochError<-sapply(1:(length(object@errors)%/%nrow(object@X)),
        function(epoch){
          sum(abs(object@errors[1:nrow(object@X)+(epoch-1)*nrow(object@X)]))
      })
      summary(epochError)
    }else{
      cat("\nStatistics for Off-line training\n")
      summary(object@errors)
    }
  } 
})
