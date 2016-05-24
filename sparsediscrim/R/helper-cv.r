#' Randomly partitions data for cross-validation.
#'
#' For a vector of training labels, we return a list of cross-validation folds,
#' where each fold has the indices of the observations to leave out in the fold.
#' In terms of classification error rate estimation, one can think of a fold as a
#' the observations to hold out as a test sample set. Either the \code{hold_out}
#' size or the number of folds, \code{num_folds}, can be specified. The number
#' of folds defaults to 10, but if the \code{hold_out} size is specified, then
#' \code{num_folds} is ignored.
#'
#' We partition the vector \code{y} based on its length, which we treat as the
#' sample size, 'n'. If an object other than a vector is used in \code{y}, its
#' length can yield unexpected results. For example, the output of
#' \code{length(diag(3))} is 9.
#'
#' @export
#' @param y a vector of class labels
#' @param num_folds the number of cross-validation folds. Ignored if
#' \code{hold_out} is not \code{NULL}. See Details.
#' @param hold_out the hold-out size for cross-validation. See Details.
#' @param seed optional random number seed for splitting the data for cross-validation
#' @return list the indices of the training and test observations for each fold.
#' @examples
#' # The following three calls to \code{cv_partition} yield the same partitions.
#' set.seed(42)
#' cv_partition(iris$Species)
#' cv_partition(iris$Species, num_folds = 10, seed = 42)
#' cv_partition(iris$Species, hold_out = 15, seed = 42)
cv_partition <- function(y, num_folds = 10, hold_out = NULL, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }
  n <- length(y)

  if (!is.null(hold_out)) {
    hold_out <- as.integer(hold_out)
    num_folds <- ceiling(n / hold_out)
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
