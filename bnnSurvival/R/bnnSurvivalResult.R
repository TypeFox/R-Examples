
setClass("bnnSurvivalResult",
  representation(
    prediction = "matrix",
    timepoints = "numeric",
    num_base_learners = "integer",
    num_features_per_base_learner = "integer",
    k = "integer",
    n_train = "integer")
)

## Constructor
bnnSurvivalResult <- function(prediction, timepoints, num_base_learners,
                              num_features_per_base_learner, k, n_train) {
  new("bnnSurvivalResult",
    prediction = prediction,
    timepoints = timepoints,
    num_base_learners = num_base_learners,
    num_features_per_base_learner = num_features_per_base_learner,
    k = k,
    n_train = n_train)
}

##' Get Predictions
##' @param object bnnSurvivalResult object to extract predictions from
##' @export
setMethod("predictions", signature("bnnSurvivalResult"),
  function(object) {
    return(object@prediction)
  }
)

##' Get timepoints
##' @param object bnnSurvivalResult object to extract timepoints from
##' @export
setMethod("timepoints", signature("bnnSurvivalResult"),
  function(object) {
    return(object@timepoints)
  }
)

##' Generic print method for bnnSurvivalResult
##' @param x bnnSurvivalResult object to print
##' @export
setMethod("print", signature("bnnSurvivalResult"),
  function(x) {
    cat("bnnSurvival results object\n\n")
    cat("Number of base learners:           ", x@num_base_learners, "\n")
    cat("Number of fatures per base learner ", x@num_features_per_base_learner, "\n")
    cat("Number of neighbors (k):           ", x@k, "\n")
    cat("Number of timepoints:              ", length(x@timepoints), "\n")
    cat("Number of training observations:   ", x@n_train, "\n")
    cat("Number of test observations:       ", nrow(x@prediction), "\n\n")
    cat("Use predictions() and timepoints() functions to access the results.\n")
  }
)

##' Generic show method for bnnSurvivalResult
##' @param object bnnSurvivalResult object to show
##' @export
setMethod("show", signature("bnnSurvivalResult"),
  function(object) {
    print(object)
  }
)
