#' PlotMF/s ANFIS' MembershipFunction domain/s
#'
#' Plot the corresponding MembershipFunctions for each/all input/s domain.
#'
#' @param object ANFIS class object.
#' @param x numeric sequence to evaluate each MembershipFunction.
#' @param input integer with the input MembershipFunctions to plot.
#' @param ... plot additional parameters.
#'
#' @return output graphics
#'
#' @include Anfis-printshow.R
#' @exportMethod plotMF
#' @docType methods
#' @name plotMF
#' @rdname ANFIS-plotMF
#' @aliases plotMF-methods
#' @seealso \code{\link{ANFIS-class}}
#' @note see full example in \code{\link{ANFIS-class}}
#' @family ANFIS
#' @author Cristobal Fresno \email{cfresno@@bdmg.com.ar}, Andrea S. Llera 
#'  \email{ALlera@@leloir.org.ar} and Elmer A. Fernandez 
#'  \email{efernandez@@bdmg.com.ar}
setGeneric(name="plotMF", def=function(object, x, input, ...){
  standardGeneric("plotMF")
})
#'
#' @name plotMF
#' @rdname ANFIS-plotMF
#' @aliases plotMF,ANFIS-method
#' @inheritParams plotMF
setMethod(f="plotMF", signature="ANFIS", 
  definition=function(object, x, input, ...){
  ##Check parameters
  if(input <= 0 | input > ncol(object@X)){
    stop("Wrong input number")
  }
  ##The domain
  plot(x, y=evaluateMF(object@premises[[input]][[1]],x), type ="l", 
    xlab=paste("Input Space #",input), ylab="Membership", ...)
  lapply(object@premises[[input]][-1],function(MF){
    lines(x,y=evaluateMF(MF,x))})
  return(invisible())
})
#'
#' @exportMethod plotMFs
#' @docType methods
#' @name plotMFs
#' @rdname ANFIS-plotMF
#' @inheritParams plotMF
#' @aliases plotMFs-methods
setGeneric(name="plotMFs", def=function(object, ...){
  standardGeneric("plotMFs")
})
#'
#' @name plotMFs
#' @rdname ANFIS-plotMF
#' @inheritParams plotMF
#' @aliases plotMFs,ANFIS-method
setMethod(f="plotMFs", signature="ANFIS", definition=function(object, ...){
  par(mfrow=c(1,ncol(object@X)))
  invisible(sapply(1:ncol(object@X),function(input){
    plotMF(object, seq(from=min(object@X[, input]), to=max(object@X[, input]),
      by=0.01), input=input, ...)
    return(NULL)
  }))
})
