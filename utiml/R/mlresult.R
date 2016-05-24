#' Convert a mlresult to a bipartition matrix
#'
#' @param mlresult The mlresult object
#' @return matrix with bipartition values
#' @export
as.bipartition <- function(mlresult) {
  utiml_ifelse(is.bipartition(mlresult),
               as.matrix(mlresult),
               attr(mlresult, "classes"))
}

#' Convert a mlresult to matrix
#'
#' @param x The mlresult object
#' @param ... ignored
#' @return matrix
#' @export
as.matrix.mlresult <- function(x, ...) {
  attr.name <- ifelse(attr(x, "type") == "bipartition", "probs", "classes")
  only.expected <- x
  attr(only.expected, attr.name) <- NULL
  attr(only.expected, "type") <- NULL
  class(only.expected) <- "matrix"
  only.expected
}

#' Convert a matrix prediction in a multi label prediction
#' @param predictions a Matrix or data.frame contained the scores/probabilities
#' values. The columns are the labels and the rows are the examples.
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#'  (Default: \code{TRUE})
#' @param ... ignored
#' @return An object of type mlresult
#' @export
#'
#' @examples
#' predictions <- matrix(runif(100), ncol = 10)
#' colnames(predictions) <- paste('label', 1:10, sep='')
#'
#' # Create a mlresult from a matrix
#' mlresult <- as.mlresult(predictions)
#' mlresult <- as.mlresult(predictions, probability = FALSE)
#' mlresult <- as.mlresult(predictions, probability = FALSE, threshold = 0.6)
#'
#' # Change the current type of a mlresult
#' mlresult <- as.mlresult(mlresult, probability = TRUE)
as.mlresult <- function(predictions, probability = TRUE, ...) {
  UseMethod("as.mlresult")
}

#' @describeIn as.mlresult Default mlresult transform method
#' @param threshold A single value between 0 and 1 or a list with threshold
#'  values contained one value per label (Default: 0.5). Only used when the
#'  predictions are not a mlresult.
#' @export
as.mlresult.default <- function (predictions, probability = TRUE, ...,
                                 threshold = 0.5) {
  predictions <- as.matrix(predictions)
  as.mlresult.mlresult(fixed_threshold(predictions, threshold), probability)
}

#' @describeIn as.mlresult change the mlresult type
#' @export
as.mlresult.mlresult <- function (predictions, probability = TRUE, ...) {
  bipartition <- as.bipartition(predictions)
  probabilities <- as.probability(predictions)
  multilabel_prediction(bipartition, probabilities, probability)
}

#' Convert a mlresult to a probability matrix
#'
#' @param mlresult The mlresult object
#' @return matrix with probabilities values
#' @export
as.probability <- function(mlresult) {
  utiml_ifelse(is.probability(mlresult),
               as.matrix(mlresult),
               attr(mlresult, "probs"))
}

#' Convert a mlresult to a ranking matrix
#'
#' @param mlresult The mlresult object
#' @param ties.method A character string specifying how ties are treated
#'  (Default: "min"). see \code{\link{rank}} to more details.
#' @param ... Others parameters passed to the \code{\link{rank}} method.
#' @return matrix with ranking values
#' @export
as.ranking <- function (mlresult, ties.method = "min", ...) {
  t(apply(1 - as.probability(mlresult), 1, rank, ties = ties.method, ...))
}

#' Test if a mlresult contains crisp values as default
#'
#' @param mlresult The mlresult object
#' @return logical value
#' @export
is.bipartition <- function(mlresult) {
  attr(mlresult, "type") == "bipartition"
}

#' Test if a mlresult contains score values as default
#'
#' @param mlresult The mlresult object
#' @return logical value
#' @export
is.probability <- function(mlresult) {
  attr(mlresult, "type") == "probability"
}

#' Create a mlresult object
#'
#' @param bipartitions The matrix of predictions (bipartition values),
#'  only 0 and 1
#' @param probabilities The matrix of probability/confidence of a prediction,
#'  between 0..1
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#'  (Default: \code{getOption("utiml.use.probs", TRUE)})
#' @return An object of type mlresult
#' @export
#' @examples
#' probs <- matrix(
#'  runif(90), ncol=3, dimnames = list(1:30, c("y1", "y2", "y3"))
#' )
#' preds <- matrix(
#'  as.numeric(probs > 0.5), ncol=3, dimnames = list(1:30, c("y1", "y2", "y3"))
#' )
#' multilabel_prediction(probs, preds)
multilabel_prediction <- function(bipartitions, probabilities,
                            probability = getOption("utiml.use.probs", TRUE)) {
  # At least one label is predict
  for (row in seq(nrow(bipartitions))) {
    bipartitions[row, probabilities[row, ] == max(probabilities[row, ])] <- 1
  }

  bipartitions <- as.matrix(bipartitions)
  probabilities <- as.matrix(probabilities)

  only.bipartitions <- bipartitions
  only.probabilities <- probabilities
  attr(probabilities, "classes") <- only.bipartitions
  attr(probabilities, "type") <- "probability"

  attr(bipartitions, "probs") <- only.probabilities
  attr(bipartitions, "type") <- "bipartition"

  class(probabilities) <- class(bipartitions) <- "mlresult"

  utiml_ifelse(probability, probabilities, bipartitions)
}

#' Print the mlresult
#' @param x The mlresult to print
#' @param ... Extra parameters for print method
#' @export
print.mlresult <- function(x, ...) {
  print(as.matrix(x), ...)
}

#' Filter a Multi-Label Result
#'
#' If column filter is performed, then the result will be a matrix. Otherwhise,
#' the result will be a mlresult.
#'
#' @param mlresult A mlresult object
#' @param rowFilter A list of rows to filter
#' @param colFilter A list of columns to filter
#' @param ... Extra parameters to be used as the filter
#' @return mlresult or matrix. If column filter is performed, then the result
#'  will be a matrix. Otherwhise, the result will be a mlresult.
#' @export
`[.mlresult` <- function (mlresult, rowFilter = T, colFilter, ...) {
  if (missing(colFilter)) {
    bipartition <- as.bipartition(mlresult)
    probability <- as.probability(mlresult)

    multilabel_prediction(bipartition[rowFilter, ],
                          probability[rowFilter, ],
                          is.probability(mlresult))
  } else {
    as.matrix(mlresult)[rowFilter, colFilter, ...]
  }
}
