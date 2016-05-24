#' Binary Relevance for multi-label Classification
#'
#' Create a Binary Relevance model for multilabel classification.
#'
#' Binary Relevance is a simple and effective transformation method to predict
#' multi-label data. This is based on the one-versus-all approach to build a
#' specific model for each label.
#'
#' @family Transformation methods
#' @param mdata A mldr dataset used to train the binary models.
#' @param base.method A string with the name of the base method. (Default:
#'  \code{options("utiml.base.method", "SVM")})
#' @param ... Others arguments passed to the base method for all subproblems
#' @param cores The number of cores to parallelize the training. Values higher
#'  than 1 require the \pkg{parallel} package. (Default:
#'  \code{options("utiml.cores", 1)})
#' @param seed An optional integer used to set the seed. This is useful when
#'  the method is run in parallel. (Default: \code{options("utiml.seed", NA)})
#' @return An object of class \code{BRmodel} containing the set of fitted
#'   models, including:
#'   \describe{
#'    \item{labels}{A vector with the label names.}
#'    \item{models}{A list of the generated models, named by the label names.}
#'   }
#' @references
#'  Boutell, M. R., Luo, J., Shen, X., & Brown, C. M. (2004). Learning
#'    multi-label scene classification. Pattern Recognition, 37(9), 1757-1771.
#' @export
#'
#' @examples
#' model <- br(toyml, "RANDOM")
#' pred <- predict(model, toyml)
#'
#' \dontrun{
#' # Use SVM as base method
#' model <- br(toyml, "SVM")
#' pred <- predict(model, toyml)
#'
#' # Change the base method and use 4 CORES
#' model <- br(toyml[1:50], 'RF', cores = 4, seed = 123)
#'
#' # Set a parameters for all subproblems
#' model <- br(toyml, 'KNN', k=5)
#' }
br <- function(mdata, base.method = getOption("utiml.base.method", "SVM"), ...,
               cores = getOption("utiml.cores", 1),
               seed = getOption("utiml.seed", NA)) {
  # Validations
  if (class(mdata) != "mldr") {
    stop("First argument must be an mldr object")
  }

  if (cores < 1) {
    stop("Cores must be a positive value")
  }

  utiml_preserve_seed()

  # BR Model class
  brmodel <- list(labels = rownames(mdata$labels), call = match.call())

  # Create models
  labels <- utiml_rename(brmodel$labels)
  brmodel$models <- utiml_lapply(labels, function (label) {
    brdata <- utiml_create_binary_data(mdata, label)
    object <- utiml_prepare_data(brdata, "mldBR", mdata$name, "br", base.method)
    utiml_create_model(object, ...)
  }, cores, seed)

  utiml_restore_seed()

  class(brmodel) <- "BRmodel"
  brmodel
}

#' Predict Method for Binary Relevance
#'
#' This function predicts values based upon a model trained by \code{\link{br}}.
#'
#' @param object Object of class '\code{BRmodel}'.
#' @param newdata An object containing the new input data. This must be a
#'  matrix, data.frame or a mldr object.
#' @param probability Logical indicating whether class probabilities should be
#'  returned. (Default: \code{getOption("utiml.use.probs", TRUE)})
#' @param ... Others arguments passed to the base method prediction for all
#'   subproblems.
#' @param cores The number of cores to parallelize the training. Values higher
#'  than 1 require the \pkg{parallel} package. (Default:
#'  \code{options("utiml.cores", 1)})
#' @param seed An optional integer used to set the seed. This is useful when
#'  the method is run in parallel. (Default: \code{options("utiml.seed", NA)})
#' @return An object of type mlresult, based on the parameter probability.
#' @seealso \code{\link[=br]{Binary Relevance (BR)}}
#' @export
#'
#' @examples
#' model <- br(toyml, "RANDOM")
#' pred <- predict(model, toyml)
#'
#' \dontrun{
#' # Predict SVM scores
#' model <- br(toyml, "SVM")
#' pred <- predict(model, toyml)
#'
#' # Predict SVM bipartitions running in 4 cores
#' pred <- predict(model, toyml, probability = FALSE, CORES = 4)
#'
#' # Passing a specif parameter for SVM predict method
#' pred <- predict(model, dataset$test, na.action = na.fail)
#' }
predict.BRmodel <- function(object, newdata,
                            probability = getOption("utiml.use.probs", TRUE),
                            ..., cores = getOption("utiml.cores", 1),
                            seed = getOption("utiml.seed", NA)) {
  # Validations
  if (class(object) != "BRmodel") {
    stop("First argument must be an BRmodel object")
  }

  if (cores < 1) {
    stop("Cores must be a positive value")
  }

  utiml_preserve_seed()

  # Create models
  newdata <- utiml_newdata(newdata)
  labels <- utiml_rename(object$labels)
  predictions <- utiml_lapply(labels, function (label) {
    utiml_predict_binary_model(object$models[[label]], newdata, ...)
  }, cores, seed)

  utiml_restore_seed()

  utiml_predict(predictions, probability)
}

#' Print BR model
#' @param x The br model
#' @param ... ignored
#' @export
print.BRmodel <- function(x, ...) {
  cat("Binary Relevance Model\n\nCall:\n")
  print(x$call)
  cat("\n", length(x$labels), "Models (labels):\n")
  print(x$labels)
}
