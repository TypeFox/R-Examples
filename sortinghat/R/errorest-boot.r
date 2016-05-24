#' Calculates the Bootstrap Error Rate for a specified classifier given a
#' data set.
#'
#' For a given data matrix and its corresponding vector of labels, we calculate
#' the bootstrap error rate for a given classifier.
#'
#' To calculate the bootstrap error rate, we sample from the data with
#' replacement to obtain a bootstrapped training data set. We then train the
#' given classifier (given in \code{train}) on the bootstrapped training data set
#' and classify the original data set given in the matrix \code{x}. Then we
#' calculate the proportion of misclassified observations, based on the true
#' labels given in \code{y}, to obtain a single bootstrap error rate. We repeat
#' this process \code{num_bootstraps} times and report the average of the
#' bootstrap error rates.
#'
#' For the given classifier, two functions must be provided 1. to train the
#' classifier and 2. to classify unlabeled observations. The training function
#' is provided as \code{train} and the classification function as
#' \code{classify}.
#'
#' We expect that the first two arguments of the \code{train} function are
#' \code{x} and \code{y}, corresponding to the data matrix and the vector of
#' their labels, respectively. Additional arguments can be passed to the
#' \code{train} function.
#'
#' We stay with the usual R convention for the \code{classify} function. We
#' expect that this function takes two arguments: 1. an \code{object} argument
#' which contains the trained classifier returned from the function specified in
#' \code{train}; and 2. a \code{newdata} argument which contains a matrix of
#' observations to be classified -- the matrix should have rows corresponding to
#' the individual observations and columns corresponding to the features
#' (covariates). For an example, see \code{\link[MASS]{lda}}.
#'
#' @export
#' @param x a matrix of n observations (rows) and p features (columns)
#' @param y a vector of n class labels
#' @param train a function that builds the classifier. (See details.)
#' @param classify a function that classifies observations from the constructed
#' classifier from \code{train}. (See details.)
#' @param num_bootstraps the number of bootstrap replications
#' @param ... additional arguments passed to the function specified in
#' \code{train}.
#' @return the bootstrapped error rate estimate
#' @examples
#' require('MASS')
#' iris_x <- data.matrix(iris[, -5])
#' iris_y <- iris[, 5]
#'
#' # Because the \code{classify} function returns multiples objects in a list,
#' # we provide a wrapper function that returns only the class labels.
#' lda_wrapper <- function(object, newdata) { predict(object, newdata)$class }
#' set.seed(42)
#' errorest_boot(x = iris_x, y = iris_y, train = MASS:::lda, classify = lda_wrapper)
#' # Output: 0.0228
errorest_boot <- function(x, y, train, classify, num_bootstraps = 50, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)
  check_out <- check_arguments(x = x, y = y, train = train, classify = classify)

  # For each bootstrap replicate, we sample with replacement from the data set,
  # construct a classifier on the bootstrap training data, classify the original
  # data in 'x', and then report the proportion of misclassified observations.
  boot_error_rates <- sapply(seq_len(num_bootstraps), function(b) {
    training <- sample(seq_along(y), replace = TRUE)
    train_out <- train(x[training, ], y[training], ...)
    classifications <- classify(object = train_out, newdata = x)
    mean(classifications != y)
  })
  
  mean(boot_error_rates)
}

