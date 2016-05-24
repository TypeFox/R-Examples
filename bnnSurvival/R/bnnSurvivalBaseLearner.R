

setClass("bnnSurvivalBaseLearner",
  representation(
    bootstrap_sample = "integer",
    feature_space = "integer")
)

## Constructor: Randomly generate bootstrap sample and feature space
bnnSurvivalBaseLearner <- function(num_samples, num_features, num_features_per_base_learner,
                                   replace, sample_fraction) {

  ## Bootstrap samples
  if (!replace & sample_fraction == 1) {
    bootstrap_sample <- 1:num_samples
  } else {
    bootstrap_sample <- sample(num_samples, num_samples * sample_fraction, replace = replace)
  }

  ## Select a subset of features if not all
  if (num_features_per_base_learner == num_features) {
    feature_space = 1:num_features
  } else {
    feature_space = sample(num_features_per_base_learner,
                           num_features_per_base_learner, replace = FALSE)
  }

  ## Create object
  new("bnnSurvivalBaseLearner",
    bootstrap_sample = bootstrap_sample,
    feature_space = feature_space)
}

##' Compute prediction for all samples.
##' @param object bnnSurvivalBaseLearner object
##' @param train_data Training data (with response)
##' @param test_data Test data (without response)
##' @param timepoints Timepoint to predict at
##' @param metric Metric used
##' @param weighting_function Weighting function used
##' @param k Number of nearest neighbors
##' @import stats 
setMethod("predict",
  signature("bnnSurvivalBaseLearner"),
  function(object, train_data, test_data, timepoints, metric, weighting_function, k) {

    ## Bootstrap sample and subsample features of training data
    train_features <- train_data[object@bootstrap_sample,
                                 object@feature_space+2, drop = FALSE]
    train_response <- train_data[object@bootstrap_sample,
                                 c(1,2), drop = FALSE]

    ## Compute distances to training obs for all test obs
    if (metric == "mahalanobis") {
      train_cov <- cov(train_features)
      
      ## Ignore all-equal features
      idx_nonzero <- rowMeans(train_cov) != 0
      
      ## Compute distances
      distances <- apply(test_data[, object@feature_space[idx_nonzero], drop = FALSE], 1,
                         mahalanobis,
                         x = train_features[, idx_nonzero, drop = FALSE],
                         cov = train_cov[idx_nonzero, idx_nonzero, drop = FALSE], 
                         tol = 1e-25)
    } else {
      stop("Currently no other distance metrics supported.")
    }

    ## Sort rows or columns, get indices
    temp <- apply(distances, 2, sort, index.return = TRUE)

    ## Compute Kaplan-Meier estimator using the k nearest neighbors for each test obs
    survival <- t(sapply(temp, function(x) {
      weighted_kaplan_meier(response = train_response[x$ix[1:k], , drop = FALSE],
                            weights = weighting_function(x$x[1:k]),
                            timepoints = timepoints)
    }))

    ## Return a matrix with predictions for all test samples and timepoints
    return(survival)
  }
)
