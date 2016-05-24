

setClass("bnnSurvivalEnsemble",
  representation(
    train_data = "matrix",
    formula = "formula",
    base_learners = "list",
    num_base_learners = "integer",
    num_features_per_base_learner = "integer",
    k = "integer",
    timepoints = "numeric",
    metric = "character",
    weighting_function = "function",
    replace = "logical", 
    sample_fraction = "numeric")
)

## Constructor
bnnSurvivalEnsemble <- function(train_data, formula, num_base_learners,
                                num_features_per_base_learner, k,
                                metric, weighting_function, 
                                replace = replace, sample_fraction = sample_fraction) {
  ## Get unique timepoints
  timepoints <- sort(unique(train_data[, 1]))

  ## Create base learners
  base_learners <- replicate(num_base_learners, bnnSurvivalBaseLearner(
                      num_samples = nrow(train_data),
                      num_features = ncol(train_data) - 2,
                      num_features_per_base_learner = num_features_per_base_learner,
                      replace = replace, sample_fraction = sample_fraction))

  new("bnnSurvivalEnsemble",
    train_data = train_data,
    formula = formula,
    base_learners = base_learners,
    num_base_learners = num_base_learners,
    num_features_per_base_learner = num_features_per_base_learner,
    k = k,
    timepoints = timepoints,
    metric = metric,
    weighting_function = weighting_function,
    replace = replace, 
    sample_fraction = sample_fraction)
}

##' Predict survival probabilities with bagged k-nearest neighbors survival prediction.
##' @param object Object of class bnnSurvivalEnsemble, created with bnnSurvival().
##' @param test_data Data set containing data to predict survival.
##' @import methods
##' @import stats
##' @importFrom parallel mclapply
##' @export
setMethod("predict", signature("bnnSurvivalEnsemble"),
  function(object, test_data) {

    ## Generate model and matrix for test data
    covar_names <- labels(terms(object@formula, data = test_data))
    test_matrix <- data.matrix(test_data[, covar_names])

    ## Check if training and test data are of same structure
    if (!all(covar_names %in% colnames(object@train_data)[c(-1, -2)])) {
      stop("Training and test data are not of same structure.")
    }

    ## Call predict on all base learners
    list_predictions <- mclapply(object@base_learners, predict, object@train_data,
                               test_matrix, object@timepoints, object@metric,
                               object@weighting_function, object@k)

    ## Aggregate predictions
    array_predictions <- simplify2array(list_predictions)
    predictions <- bnnSurvivalPredictions(array_predictions, object@timepoints, object@num_base_learners,
                                          object@num_features_per_base_learner, object@k,
                                          nrow(object@train_data))
    result <- aggregate(predictions)

    ## Return result
    return(result)
  }
)

##' Generic print method for bnnSurvivalEnsemble
##' @param x bnnSurvivalEnsemble object to print
##' @import methods
##' @export
setMethod("print", signature("bnnSurvivalEnsemble"),
  function(x) {
    cat("bnnSurvival ensemble object\n\n")
    cat("Formula:                           ", deparse(x@formula), "\n")
    cat("Number of base learners:           ", x@num_base_learners, "\n")
    cat("Number of fatures per base learner ", x@num_features_per_base_learner, "\n")
    cat("Number of neighbors (k):           ", x@k, "\n")
    cat("Number of timepoints:              ", length(x@timepoints), "\n")
    cat("Number of training observations:   ", nrow(x@train_data), "\n")
    cat("Used metric:                       ", x@metric, "\n")
    cat("Weighting function:                ", deparse(x@weighting_function), "\n")
    cat("Sample with replacement:           ", x@replace, "\n")
    cat("Sample fraction:                   ", x@sample_fraction, "\n\n")
    cat("Use predict() method to predict surival probabilities for new data.\n")
  }
)

##' Generic show method for bnnSurvivalEnsemble
##' @param object bnnSurvivalEnsemble object to show
##' @import methods
##' @export
setMethod("show", signature("bnnSurvivalEnsemble"),
  function(object) {
    print(object)
  }
)
