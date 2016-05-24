
##' Get Predictions
##' @param object Object to extract predictions from
##' @param ... further arguments passed to or from other methods.
setGeneric("predictions", function(object, ...) standardGeneric("predictions"))

##' Get Timepoints
##' @param object Object to extract timepoints from
##' @param ... further arguments passed to or from other methods.
setGeneric("timepoints", function(object, ...) standardGeneric("timepoints"))

##' Function to extract survival probability predictions from bnnSurvivalEnsemble. Use with \code{pec} package.
##' @param object bnnSurvivalEnsemble object.
##' @param newdata Data used for prediction.
##' @param times Not used.
##' @param ... Not used.
##' @return survival probability predictions
##' @importFrom pec predictSurvProb
##' @export
predictSurvProb.bnnSurvivalEnsemble <- function(object, newdata, times, ...) {
  pred <- predict(object, newdata)
  pos <- sapply(times, function(x) {
    max(c(0, which(x >= timepoints(pred))))
  })    
  p <- cbind(1, predictions(pred))[, pos+1, drop = FALSE]
  return(p)
}