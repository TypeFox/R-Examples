#' Partitions data for cross-validation.
#'
#' For a vector of training labels, we return a list of cross-validation folds,
#' where each fold has the indices of the observations to leave out in the fold.
#' In terms of classification error rate estimation, one can think of a fold as
#' a the observations to hold out as a test sample set.
#'
#' Either the \code{hold_out} size or \code{num_folds} can be specified. The
#' number of folds defaults to 10, but if the \code{hold_out} size is specified,
#' then \code{num_folds} is ignored.
#'
#' We partition the vector \code{y} based on its length, which we treat as the
#' sample size, \code{n}. If an object other than a vector is used in \code{y},
#' its length can yield unexpected results. For example, the output of
#' \code{length(diag(3))} is 9.
#'
#' @export
#' @param y a vector of class labels to partition
#' @param num_folds the number of cross-validation folds. Ignored if
#' \code{hold_out} is not \code{NULL}. See Details.
#' @param hold_out the hold-out size for cross-validation. See Details.
#' @param seed optional random number seed for splitting the data for cross-validation
#' @return list the indices of the training and test observations for each fold.
#' @examples
#' library(MASS)
#' # The following three calls to \code{cv_partition} yield the same partitions.
#' set.seed(42)
#' cv_partition(iris$Species)
#' cv_partition(iris$Species, num_folds = 10, seed = 42)
#' cv_partition(iris$Species, hold_out = 15, seed = 42)
cv_partition <- function(y, num_folds = 10, hold_out = NULL, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n <- length(y)

  if (!is.null(hold_out)) {
    hold_out <- as.integer(hold_out)
    num_folds <- ceiling(n / hold_out)
  } else if (num_folds > n) {
    warning("The number of folds exceeds length(y). Setting 'num_folds' to 'n'...")
    num_folds <- n
  }

  folds <- split(sample(seq_len(n), n), gl(n = num_folds, k = 1, length = n))
  folds <- lapply(folds, function(fold) {
    list(
      training = which(!seq_along(y) %in% fold),
      test = fold
    )
  })
  names(folds) <- paste0("Fold", names(folds))

  folds
}

#' Helper function that partitions a data set into training and test data sets.
#'
#' The function randomly partitions a data set into training and test data sets
#' with a specified percentage of observations assigned to the training data set.
#' The user can optionally preserve the proportions of the original data set.
#'
#' A named list is returned with the training and test data sets.
#'
#' @export
#' @param x a matrix of n observations (rows) and p features (columns)
#' @param y a vector of n class labels
#' @param split_pct the percentage of observations that will be randomly assigned
#' to the training data set. The remainder of the observations will be assigned
#' to the test data set.
#' @param preserve_proportions logical value. If \code{TRUE}, the training and
#' test data sets will be constructed so that the original proportions are
#' preserved.
#' @return named list containing the training and test data sets:
#' \itemize{
#'   \item \code{train_x}: matrix of the training observations
#'   \item \code{train_y}: vector of the training labels (coerced to factors).
#'   \item \code{test_x}: matrix of the test observations
#'   \item \code{test_y}: vector of the test labels (coerced to factors).
#' }
#' @examples
#' require('MASS')
#' x <- iris[, -5]
#' y <- iris[, 5]
#' set.seed(42)
#' data <- partition_data(x = x, y = y)
#' table(data$train_y)
#' table(data$test_y)
#'
#' data <- partition_data(x = x, y = y, preserve_proportions = TRUE)
#' table(data$train_y)
#' table(data$test_y)
partition_data <- function(x, y, split_pct = 2/3, preserve_proportions = FALSE) {
  x <- as.matrix(x)
  y <- as.factor(y)

  if (preserve_proportions) {
    train <- tapply(seq_along(y), y, function(i) {
      sample(i, size = length(i) * split_pct)
    })
    train <- do.call(c, train)
  } else {
    train <- sample(seq_along(y), size = length(y) * split_pct)
  }
  list(train_x = x[train, ], train_y = y[train], test_x = x[-train, ],
       test_y = y[-train])
}

