

setClass("bnnSurvivalPredictions",
  representation(
    predictions = "array", 
    timepoints = "numeric",
    num_base_learners = "integer",  
    num_features_per_base_learner = "integer",
    k = "integer", 
    n_train = "integer")
)

## Constructor
bnnSurvivalPredictions <- function(predictions, timepoints, num_base_learners, 
                                   num_features_per_base_learner, k, n_train) {
  new("bnnSurvivalPredictions", 
      predictions = predictions,
      timepoints = timepoints, 
      num_base_learners = num_base_learners,
      num_features_per_base_learner = num_features_per_base_learner,
      k = k,
      n_train = n_train)
}

setMethod("aggregate",
  signature("bnnSurvivalPredictions"),
    function(x) {
      ## Aggregate all predictions
      aggregated_predictions <- apply(x@predictions, c(1,2), mean, 
                                      na.rm = TRUE)
      
      ## Create and return Result object
      result <- bnnSurvivalResult(aggregated_predictions, x@timepoints, x@num_base_learners, 
                                  x@num_features_per_base_learner, x@k, x@n_train)
      return(result)   
    }
)