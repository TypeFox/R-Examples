#' Compute the multi-label ensemble predictions based on some vote schema
#'
#' @param predictions A list of multi-label predictions (mlresult).
#' @param vote.schema Define the way that ensemble must compute the predictions.
#'  The default valid options are:
#'  \describe{
#'    \item{'avg'}{Compute the mean of probabilities and the bipartitions}
#'    \item{'maj'}{Compute the majority of votes}
#'    \item{'max'}{Compute the higher probability for each instance/label}
#'    \item{'min'}{Compute the lower probability for each instance/label}
#'  }. (Default: 'maj')
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#' @return A mlresult with computed predictions.
#' @note You can create your own vote schema, just create a method that receive
#'  two matrix (bipartitions and probabilities) and return a list with the
#'  final bipartions and probabilities.
#'
#'  Remember that this method will compute the ensemble votes for each label.
#'  Thus the bipartition and probability matrix passed as argument for this
#'  method is related with the bipartitions and probabilities for a single
#'  label.
#' @export
#'
#' @examples
#' \dontrun{
#' model <- br(toyml, "KNN")
#' predictions <- list(
#'  predict(model, toyml[1:10], k=1),
#'  predict(model, toyml[1:10], k=3),
#'  predict(model, toyml[1:10], k=5)
#' )
#'
#' result <- compute_multilabel_predictions(predictions, "maj")
#'
#' ## Random choice
#' random_choice <- function (bipartition, probability) {
#'  cols <- sample(seq(ncol(bipartition)), nrow(bipartition), replace = TRUE)
#'  list(
#'    bipartition = bipartition[cbind(seq(nrow(bipartition)), cols)],
#'    probability = probability[cbind(seq(nrow(probability)), cols)]
#'  )
#' }
#' result <- compute_multilabel_predictions(predictions, "random_choice")
#' }
compute_multilabel_predictions <- function(predictions,
                                           vote.schema = "maj",
                                probability = getOption("utiml.use.probs", TRUE)
                                           ) {

  if (length(unique(lapply(predictions, dimnames))) > 1) {
    stop("The predictions must be the same dimensions and names.")
  }

  utiml_ensemble_check_voteschema(vote.schema, FALSE)
  vote.method <- utiml_ensemble_method(vote.schema)

  probs <- lapply(predictions, as.probability)
  preds <- lapply(predictions, as.bipartition)

  examples <- rownames(probs[[1]])
  labels <- utiml_rename(colnames(probs[[1]]))
  final.prediction <- lapply(labels, function (label) {
    lpreds <- do.call(cbind,
                      lapply(preds, function (prediction) prediction[, label]))
    lprobs <- do.call(cbind,
                      lapply(probs, function (prediction) prediction[, label]))

    utiml_compute_ensemble(lpreds, lprobs, vote.method, examples)
  })

  utiml_predict(final.prediction, probability)
}

# Internal methods -------------------------------------------------------------
#' Average vote combination for a single-label prediction
#'
#' Compute the prediction for a single-label using the average votes schema.
#' The probabilities result is computed using the averaged values.
#'
#' @param bipartition A matrix with all bipartition predictions for a single
#'  label. The column are the predictions and the rows the examples.
#' @param probability A matrix with all probability predictions for a single
#'  label The column are the predictions and the rows the examples.
#' @return A list with two values "bipartition" and "probability".
utiml_ensemble_average <- function (bipartition, probability) {
  list(
    bipartition = as.numeric(rowMeans(bipartition) >= 0.5),
    probability = rowMeans(probability)
  )
}

#' Verify if a schema vote name is valid
#'
#' @param vote.schema The name of schema vote
#' @param accept.null Logical value determine if the vote.schema = NULL is
#'  also valid. (Default: TRUE)
#' @return TRUE or throw an error message otherwise
utiml_ensemble_check_voteschema <- function (vote.schema, accept.null = TRUE) {
  if (is.null(vote.schema)) {
    if (!accept.null) {
      stop("The enseble vote schema can not be NULL")
    }
  }
  else if (!vote.schema %in% c("avg", "maj", "max", "min")) {
    if (!exists(vote.schema, mode = "function")) {
      stop(paste("The compute ensemble method '", vote.schema,
                 "' is not a valid function", sep=''))
    }
  }
  invisible(TRUE)
}

#' Majority vote combination for single-label prediction
#'
#' Compute the single-label prediction using the majority votes schema.
#' The probabilities result is computed using only the majority instances.
#' In others words, if a example is predicted as posivite, only the positive
#' confidences are used to compute the averaged value.
#'
#' @param bipartition A matrix with all bipartition predictions for a single
#'  label. The column are the predictions and the rows the examples.
#' @param probability A matrix with all probability predictions for a single
#'  label The column are the predictions and the rows the examples.
#' @return A list with two values "bipartition" and "probability".
utiml_ensemble_majority <- function (bipartition, probability) {
  # Compute the votes
  votes <- rowMeans(bipartition)
  ties <- votes == 0.5
  votes[ties] <- rowMeans(as.matrix(probability[ties,]))
  votes <- as.numeric(votes >= 0.5)

  probs <- votes
  positive <- votes == 1

  # Compute the positive probabilities
  probs[positive] <- unlist(lapply(which(positive), function(row) {
    mean(probability[row, bipartition[row, ] == 1])
  }))

  # Compute the negative p
  probs[!positive] <- unlist(lapply(which(!positive), function(row) {
    mean(probability[row, bipartition[row, ] == 0])
  }))

  list(bipartition = votes, probability = probs)
}

#' Maximum vote combination for single-label prediction
#'
#' Compute the single-label prediction using the maximum votes schema. The
#' probabilities result is computed using the maximum value.
#'
#' @param bipartition A matrix with all bipartition predictions for a single
#'  label. The column are the predictions and the rows the examples.
#' @param probability A matrix with all probability predictions for a single
#'  label The column are the predictions and the rows the examples.
#' @return A list with two values "bipartition" and "probability".
utiml_ensemble_maximum <- function (bipartition, probability) {
  list(
    bipartition = apply(bipartition, 1, max),
    probability = apply(probability, 1, max)
  )
}

#' Define the method name related with the vote schema
#'
#' @param vote.schema Define the way that ensemble must compute the predictions.
#' @return The method name that will compute the votes
utiml_ensemble_method <- function(vote.schema) {
  votes <- c(
    avg  = "utiml_ensemble_average",
    maj  = "utiml_ensemble_majority",
    max  = "utiml_ensemble_maximum",
    min  = "utiml_ensemble_minimum"
  )

  as.character(ifelse(is.na(votes[vote.schema]),
                      vote.schema, votes[vote.schema]))
}

#' Minimum vote combination for single-label prediction
#'
#' Compute the single-label prediction using the minimum votes schema. The
#' probabilities result is computed using the minimum value.
#'
#' @param bipartition A matrix with all bipartition predictions for a single
#'  label. The column are the predictions and the rows the examples.
#' @param probability A matrix with all probability predictions for a single
#'  label The column are the predictions and the rows the examples.
#' @return A list with two values "bipartition" and "probability".
utiml_ensemble_minimum <- function (bipartition, probability) {
  list(
    bipartition = apply(bipartition, 1, min),
    probability = apply(probability, 1, min)
  )
}

#' @describeIn compute_multilabel_predictions Internal version
utiml_predict_ensemble <- function(predictions, vote.schema,
                                            probability) {
  if (is.null(vote.schema)) {
    return(predictions)
  } else {
    compute_multilabel_predictions(predictions, vote.schema, probability)
  }
}

#' Compute binary predictions
#'
#' @param bipartitions A matrix with bipartitions values.
#' @param probabilities A matrix with probabilities values.
#' @param vote.methods The vote schema method.
#' @param rnames The row names.
#'
#' @return A binary.prediction object
utiml_compute_ensemble <- function (bipartitions, probabilities,
                                    vote.methods, rnames) {
  result <- do.call(vote.methods, list(bipartitions, probabilities))
  names(result$bipartition) <- names(result$probability) <- rnames
  utiml_binary_prediction(result$bipartition, result$probability)
}

#' Predict binary predictions
#'
#' Is very simillar from utiml_compute_ensemble but differs from arguments
#'
#' @param predictions A list of binary predictions.
#' @param vote.schema The name of vote schema.
#'
#' @return A binary.prediction object
utiml_predict_binary_ensemble <- function(predictions, vote.schema) {
  lpreds <- do.call(cbind,
                    lapply(predictions, function (pred) pred$bipartition))
  lprobs <- do.call(cbind,
                    lapply(predictions, function (pred) pred$probability))

  utiml_compute_ensemble(lpreds, lprobs,
                         utiml_ensemble_method(vote.schema),
                         rownames(lpreds))
}
