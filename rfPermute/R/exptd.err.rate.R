#' @title Expected Error Rate
#' @description Calculate expected OOB error rates for randomForest 
#'   classification model based on random assignment and class sizes (prior).
#' 
#' @param rf an object inheriting from \code{link{randomForest}}.
#' 
#' @return a vector of expected error rates for each class.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
exptd.err.rate <- function(rf) {
  if(!inherits(rf, "randomForest")) {
    stop("'rf' is not a randomForest or rfPermute object.")
  }
  if(rf$type == "regression") stop("'rf' must be of a classification model")
  y.freq <- table(rf$y)
  exp.err <- 1 - y.freq / sum(y.freq)
  oob.err <- sum(exp.err * y.freq) / sum(y.freq)
  c(OOB = oob.err, exp.err)
}
