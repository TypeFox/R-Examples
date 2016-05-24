######################################################################################################################

# Function: CreateEvaluationStructure.
# Argument: Evaluation model.
# Description: This function is based on the old evaluation_model_extract function. It performs error checks in the evaluation model
# and creates an "evaluation structure", which is an internal representation of the original evaluation model used by all other Mediana functions.

CreateEvaluationStructure = function(evaluation.model) {

  # TO DO Make sure that all criteria IDs are different

  # General set of evaluation model parameters
  general = evaluation.model$general

  if (is.null(evaluation.model$criteria))
    stop("Evaluation model: At least one criterion must be specified.")

  # Extract criterion-specific parameters

  # Number of criteria
  n.criteria = length(evaluation.model$criteria)

  # List of criteria (id, method, test list, statistic list, parameters, result label list)
  criterion = list()

  for (i in 1:n.criteria) {
    # Metric IDs
    if (is.null(evaluation.model$criteria[[i]]$id))
      stop("Evaluation model: IDs must be specified for all criteria.") else id = evaluation.model$criteria[[i]]$id
      # Criteria
      if (is.null(evaluation.model$criteria[[i]]$method)) {
        stop("Evaluation model: Criterion method must be specified for all criteria.")
      } else if (!exists(evaluation.model$criteria[[i]]$method)) {
        stop(paste0("Evaluation model: Criterion function '", evaluation.model$criteria[[i]]$method, "' does not exist."))
      } else if (!is.function(get(as.character(evaluation.model$criteria[[i]]$method), mode = "any"))) {
        stop(paste0("Evaluation model: Criterion function '", evaluation.model$criteria[[i]]$method, "' does not exist."))
      } else {
        method = evaluation.model$criteria[[i]]$method
      }

      # Tests and statistics
      if (is.null(evaluation.model$criteria[[i]]$tests) & is.null(evaluation.model$criteria[[i]]$statistics))
        stop("Evaluation model: Tests or statistics must be specified for all criteria.")

      if (!is.null(evaluation.model$criteria[[i]]$tests)) {
        tests = evaluation.model$criteria[[i]]$tests
      } else {
        tests = NULL
      }

      if (!is.null(evaluation.model$criteria[[i]]$statistics)) {
        statistics = evaluation.model$criteria[[i]]$statistics
      } else {
        statistics = NULL
      }

      # Parameters (optional)
      if (is.null(evaluation.model$criteria[[i]]$par)) {
        par = NA
      } else {
        par = evaluation.model$criteria[[i]]$par
      }

      # Result labels
      if (is.null(evaluation.model$criteria[[i]]$labels)) {
        stop(paste0("Evaluation model: Label must be specified for the criterion ",evaluation.model$criteria[[i]]$id,"."))
      } else {
        labels = evaluation.model$criteria[[i]]$labels
      }

      criterion[[i]] = list(id = id, method = method, tests = tests, statistics = statistics, par = par, labels = labels)
  }

  # Create the evaluation structure
  evaluation.structure = list(description = "evaluation.structure",
                              criterion = criterion,
                              general = general)
  return(evaluation.structure)

}
# End of CreateEvaluationStructure