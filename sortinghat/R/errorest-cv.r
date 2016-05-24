#' Calculates the Cross-Validation Error Rate for a specified classifier given a
#' data set.
#'
#' For a given data matrix and its corresponding vector of labels, we calculate
#' the cross-validation (CV) error rate for a given classifier.
#'
#' To calculate the CV error rate, we partition the data set into 'folds'. For
#' each fold, we consider the observations within the fold as a test data set,
#' while the remaining observations are considered as a training data set. We
#' then calculate the number of misclassified observations within the fold. The
#' CV error rate is the proportion of misclassified observations across all
#' folds.
#'
#' Rather than partitioning the observations into folds, an alternative
#' convention is to specify the 'hold-out' size for each test data set. Note that
#' this convention is equivalent to the notion of folds. We allow the user to
#' specify either option with the \code{hold_out} and \code{num_folds} arguments.
#' The \code{num_folds} argument is the default option but is ignored if the
#' \code{hold_out} argument is specified (i.e. is not \code{NULL}).
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
#' @param num_folds the number of cross-validation folds. Ignored if
#' \code{hold_out} is not \code{NULL}. See Details.
#' @param hold_out the hold-out size for cross-validation. See Details.
#' @param ... additional arguments passed to the function specified in
#' \code{train}.
#' @return the calculated CV error-rate estimate
#' @examples
#' require('MASS')
#' iris_x <- data.matrix(iris[, -5])
#' iris_y <- iris[, 5]
#'
#' # Because the \code{classify} function returns multiples objects in a list,
#' # we provide a wrapper function that returns only the class labels.
#' lda_wrapper <- function(object, newdata) { predict(object, newdata)$class }
#' 
#' set.seed(42)
#' errorest_cv(x = iris_x, y = iris_y, train = MASS:::lda, classify = lda_wrapper)
#' # Output: 0.02666667
errorest_cv <- function(x, y, train, classify, num_folds = 10, hold_out = NULL, ...) {
  x <- as.matrix(x)
  y <- as.factor(y)
  check_out <- check_arguments(x = x, y = y, train = train, classify = classify)

  list_partitions <- cv_partition(y, num_folds = num_folds, hold_out = hold_out)

  # For each fold, we calculate the number of misclassified observations.
  num_misclassified <- sapply(list_partitions, function(fold) {
    train_out <- with(fold, train(x[training, ], y[training], ...))
    classifications <- with(fold, classify(object = train_out,
                                          newdata = x[-training, ]))
    with(fold, sum(classifications != y[-training]))
  })
  # To compute the CV error rate, we sum up the total number of misclassified
  # observations and then divide through by the sample size (the length of 'y').
  sum(num_misclassified) / length(y)
}

