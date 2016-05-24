.empRocInfo <- function(scores, classes) {
  # This software comes with absolutely no warranty. Use at your own risk.
  #
  # Provides information related to the ROC given probability and class label vectors. 
  # This function is not to be called directly in a normal use case. 
  # Instead, the other functions in this package call this function when necessary.
  #
  #
  # Arguments:
  #   scores: A vector of probability scores.
  #   classes: A vector of true class labels.
  # Value:
  #   A RocInfo object with six components:
  #     n0: Number of positive observations.
  #     n1: Number of negative observations.
  #     pi0: Prior probability of positive observation.
  #     pi1: Prior probability of negative observation.
  #     F0: Convex hull of ROC y values.
  #     F1: Convex hull of ROC x values.
  if (length(scores) != length(classes)) {
    stop('Length of scores and classes vectors is not equal')
  }
  prediction <- prediction(scores, classes)
  perf <- performance(prediction, "rch")
  n0 <- prediction@n.pos[[1]]
  n1 <- prediction@n.neg[[1]]
  pi0 <- n0 / (n0 + n1)
  pi1 <- n1 / (n0 + n1)
  F0 <- perf@y.values[[1]]
  F1 <- perf@x.values[[1]]
  list(n0=n0, n1=n1, pi0=pi0, pi1=pi1, F0=F0, F1=F1)
}
