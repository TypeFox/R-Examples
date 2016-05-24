#' Calculates the Apparent Error Rate for a specified classifier given a data
#' set.
#'
#' For a given data matrix and its corresponding vector of labels, we calculate
#' the apparent error rate (AER) for a given classifier.
#'
#' The AER simply uses the data set as both the training and test data sets. The
#' AER is well known to be biased downward in that it underestimates the true
#' error rate of the classifier.
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
#' @param classify a function that classifies observations from the classifier
#' constructed by \code{train}. (See details.)
#' @param ... additional arguments passed to the function specified in
#' \code{train}.
#' @return object of class \code{errorest}. The object is a named \code{list}
#' that contains the following elements:
#' @return the calculated apparent error-rate estimate
#' @examples
#' require('MASS')
#' iris_x <- data.matrix(iris[, -5])
#' iris_y <- iris[, 5]
#'
#' # Because the \code{classify} function returns multiples objects in a list,
#' # we provide a wrapper function that returns only the class labels.
#' lda_wrapper <- function(object, newdata) { predict(object, newdata)$class }
#' errorest_apparent(x = iris_x, y = iris_y, train = MASS:::lda, classify = lda_wrapper)
#' # Output: 0.02
#' 
#' # The following code is equivalent for this example:
#' lda_out <- MASS:::lda(x = iris_x, grouping = iris_y)
#' lda_classifications <- predict(lda_out, newdata = iris_x)$class
#' mean(lda_classifications != iris_y)
#' # Output: 0.02
errorest_apparent <- function(x, y, train, classify, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)
  check_out <- check_arguments(x = x, y = y, train = train, classify = classify)

  train_out <- train(x, y, ...)
  classifications <- classify(object = train_out, newdata = x)
  mean(y != classifications)
}
