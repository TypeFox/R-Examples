
## NAMESPACE

#' @importFrom randomForest randomForest
NULL

#' @importFrom methods setClass setGeneric setMethod  extends getClass is new prototype signature slot
NULL

#' @import caret
NULL

#' @importFrom stats cor lowess predict quantile rbinom sd
NULL

#' @importFrom utils tail head
NULL

## FUNCTIONS

# Mode for VOTE

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]}

# Test for even/odd

is.odd <- function(x) x %% 2 != 0

# Extract validation results

extract <- function(x){
  row <- c(variance=x@variance, finite=x@finite, completeobs=x@completeobs, classbalance=x@classbalance, ntopratiotwoplus=x@ntopratiotwoplus, mindimensions=x@mindimensions)
}



#globalVariables(c("result","combinationevaluation", "predictor", "skewness"))

