#' Classifier Chains for multi-label Classification
#'
#' Create a Classifier Chains model for multilabel classification.
#'
#' Classifier Chains is a Binary Relevance transformation method based to
#' predict multi-label data. This is based on the one-versus-all approach to
#' build a specific model for each label. It is different from BR method due the
#' strategy of extended the attribute space with the 0/1 label relevances of all
#' previous classifiers, forming a classifier chain.
#'
#' @family Transformation methods
#' @param mdata A mldr dataset used to train the binary models.
#' @param base.method A string with the name of the base method. (Default:
#'  \code{options("utiml.base.method", "SVM")})
#' @param chain A vector with the label names to define the chain order. If
#'   empty the chain is the default label sequence of the dataset. (Default:
#'   \code{NA})
#' @param ... Others arguments passed to the base method for all subproblems.
#' @param cores The number of cores to parallelize the training. Values higher
#'  than 1 require the \pkg{parallel} package. (Default:
#'  \code{options("utiml.cores", 1)})
#' @param seed An optional integer used to set the seed. This is useful when
#'  the method is run in parallel. (Default: \code{options("utiml.seed", NA)})
#' @return An object of class \code{CCmodel} containing the set of fitted
#'   models, including: \describe{
#'   \item{chain}{A vector with the chain order.}
#'   \item{labels}{A vector with the label names in expected order.}
#'   \item{models}{A list of models named by the label names.}
#' }
#' @references
#'  Read, J., Pfahringer, B., Holmes, G., & Frank, E. (2011). Classifier chains
#'    for multi-label classification. Machine Learning, 85(3), 333-359.
#'
#'  Read, J., Pfahringer, B., Holmes, G., & Frank, E. (2009). Classifier Chains
#'    for Multi-label Classification. Machine Learning and Knowledge Discovery
#'    in Databases, Lecture Notes in Computer Science, 5782, 254-269.
#' @export
#'
#' @examples
#' model <- cc(toyml, "RANDOM")
#' pred <- predict(model, toyml)
#'
#' \dontrun{
#' # Use a specific chain with J48 classifier
#' mychain <- sample(rownames(toyml$labels))
#' model <- cc(toyml, 'J48', mychain)
#'
#' # Set a specific parameter
#' model <- cc(toyml, 'KNN', k=5)
#'
#' #Run with multiple-cores
#' model <- cc(toyml, 'RF', cores = 5, seed = 123)
#' }
cc <- function(mdata, base.method = getOption("utiml.base.method", "SVM"),
               chain = NA, ..., cores = getOption("utiml.cores", 1),
               seed = getOption("utiml.seed", NA)) {
  # Validations
  if (class(mdata) != "mldr") {
    stop("First argument must be an mldr object")
  }

  labels <- rownames(mdata$labels)
  chain <- utiml_ifelse(anyNA(chain), labels, chain)
  if (!utiml_is_equal_sets(chain, labels)) {
    stop("Invalid chain (all labels must be on the chain)")
  }

  utiml_preserve_seed()

  # CC Model class
  ccmodel <- list(labels = labels, chain = chain, call = match.call())

  # Create models
  basedata <- mdata$dataset[mdata$attributesIndexes]
  labeldata <- mdata$dataset[mdata$labels$index][chain]
  chain.order <- utiml_rename(seq(mdata$measures$num.labels), chain)
  ccmodel$models <- utiml_lapply(chain.order, function(lidx) {
    data <- cbind(basedata, labeldata[seq(lidx)])
    dataset <- utiml_prepare_data(data, "mldCC", mdata$name, "cc", base.method,
                               chain.order = lidx)
    utiml_create_model(dataset, ...)
  }, cores, seed)

  utiml_restore_seed()

  class(ccmodel) <- "CCmodel"
  ccmodel
}

#' Predict Method for Classifier Chains
#'
#' This function predicts values based upon a model trained by \code{cc}.
#'
#' @param object Object of class '\code{CCmodel}'.
#' @param newdata An object containing the new input data. This must be a
#'  matrix, data.frame or a mldr object.
#' @param probability Logical indicating whether class probabilities should be
#'  returned. (Default: \code{getOption("utiml.use.probs", TRUE)})
#' @param ... Others arguments passed to the base method prediction for all
#'   subproblems.
#' @param cores Ignored because this method does not support multi-core.
#' @param seed An optional integer used to set the seed.
#'  (Default: \code{options("utiml.seed", NA)})
#' @return An object of type mlresult, based on the parameter probability.
#' @seealso \code{\link[=cc]{Classifier Chains (CC)}}
#' @note The Classifier Chains prediction can not be parellelized.
#' @export
#'
#' @examples
#' model <- cc(toyml, "RANDOM")
#' pred <- predict(model, toyml)
#'
#' \dontrun{
#' # Predict SVM bipartitions
#' pred <- predict(model, toyml, prob = FALSE)
#'
#' # Passing a specif parameter for SVM predict method
#' pred <- predict(model, toyml, na.action = na.fail)
#' }
predict.CCmodel <- function(object, newdata,
                            probability = getOption("utiml.use.probs", TRUE),
                            ..., cores = NULL,
                            seed = getOption("utiml.seed", NA)) {
  # Validations
  if (class(object) != "CCmodel") {
    stop("First argument must be an CCmodel object")
  }

  utiml_preserve_seed()
  if (!anyNA(seed)) {
    set.seed(seed)
  }

  newdata <- list(utiml_newdata(newdata))
  predictions <- list()
  for (label in object$chain) {
    predictions[[label]] <- utiml_predict_binary_model(object$models[[label]],
                                                       do.call(cbind, newdata),
                                                       ...)
    newdata[[label]] <- predictions[[label]]$bipartition
  }

  utiml_restore_seed()
  utiml_predict(predictions[object$labels], probability)
}

#' Print CC model
#' @param x The cc model
#' @param ... ignored
#' @export
print.CCmodel <- function(x, ...) {
  cat("Classifier Chains Model\n\nCall:\n")
  print(x$call)
  cat("\nChain: (", length(x$chain), "labels )\n")
  print(paste(x$chain, collapse =' -> '))
}
