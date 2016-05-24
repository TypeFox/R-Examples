#' Methods for function \code{addWeights}
#'
#' allows to modify sampling weights of an \code{\linkS4class{dataObj}} or
#' \code{\linkS4class{simPopObj}}-object. As input the output of
#' \code{\link{calibSample}} must be used.
#' @rdname addWeights
#' @aliases addWeights
#' @param object an object of class \code{\linkS4class{dataObj}} or \code{\linkS4class{simPopObj}}.
#' @param value a numeric vector of suitable length
#' @export
#' @examples
#' data(eusilcS)
#' data(totalsRG)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' \dontrun{
#' ## approx. 20 seconds ...
#' addWeights(inp) <- calibSample(inp, totalsRG)
#' }
setGeneric("addWeights<-", function(object, value) standardGeneric("addWeights<-"))

#' @aliases addWeights<-,dataObj-method
#' @rdname addWeights
#' @export
setReplaceMethod("addWeights", signature = c("dataObj"), definition=function(object, value) {
  if ( length(value) != 2 ) {
    stop("The provided input must be the output of 'calibWeights()'\n")
  }
  if ( sum(names(value) != c("g_weights", "final_weights")) > 0 ) {
    stop("The provided input is not valid. It must be an output of 'calibWeights()'!\n")
  }
  if ( length(value[[1]]) != nrow(object@data) ) {
    stop("The dimensions do not match. Please check your input!\n")
  }
  object@data[[object@weight]] <- value$final_weights
  validObject(object)
  invisible(object)
})

#' @aliases addWeights<-,simPopObj-method
#' @rdname addWeights
#' @export
setReplaceMethod("addWeights", signature = c("simPopObj"), definition=function(object, value) {
  if ( is.null(object@sample) ) {
    stop("No sample information is provided in the input object!\n")
  }
  dat <- object@sample
  addWeights(dat) <- value
  object@sample <- dat
  validObject(object)
  invisible(object)
})

