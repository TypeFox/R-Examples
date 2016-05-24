#' Ensemble of Classifier Chains for multi-label Classification
#'
#' Create an Ensemble of Classifier Chains model for multilabel classification.
#'
#' This model is composed by a set of Classifier Chains models. Classifier
#' Chains is a Binary Relevance transformation method based to predict
#' multi-label data. It is different from BR method due the strategy of extended
#' the attribute space with the 0/1 label relevances of all previous
#' classifiers, forming a classifier chain.
#'
#' @family Transformation methods
#' @family Ensemble methods
#' @param mdata A mldr dataset used to train the binary models.
#' @param base.method A string with the name of the base method. (Default:
#'  \code{options("utiml.base.method", "SVM")})
#' @param m The number of Classifier Chains models used in the ensemble.
#' @param subsample A value between 0.1 and 1 to determine the percentage of
#'    training instances that must be used for each classifier. (Default: 0.75)
#' @param attr.space A value between 0.1 and 1 to determine the percentage of
#'    attributes that must be used for each classifier. (Default: 0.50)
#' @param ... Others arguments passed to the base method for all subproblems.
#' @param cores The number of cores to parallelize the training. Values higher
#'  than 1 require the \pkg{parallel} package. (Default:
#'  \code{options("utiml.cores", 1)})
#' @param seed An optional integer used to set the seed. This is useful when
#'  the method is run in parallel. (Default: \code{options("utiml.seed", NA)})
#' @return An object of class \code{ECCmodel} containing the set of fitted
#'   CC models, including:
#' \describe{
#'   \item{rounds}{The number of interations}
#'   \item{models}{A list of BR models.}
#'   \item{nrow}{The number of instances used in each training dataset}
#'   \item{ncol}{The number of attributes used in each training dataset}
#' }
#' @references
#'    Read, J., Pfahringer, B., Holmes, G., & Frank, E. (2011). Classifier
#'    chains for multi-label classification. Machine Learning, 85(3), 333-359.
#'
#'    Read, J., Pfahringer, B., Holmes, G., & Frank, E. (2009).
#'    Classifier Chains for Multi-label Classification. Machine Learning and
#'    Knowledge Discovery in Databases, Lecture Notes in Computer Science,
#'    5782, 254-269.
#' @note If you want to reproduce the same classification and obtain the same
#'  result will be necessary set a flag utiml.mc.set.seed to FALSE.
#' @export
#'
#' @examples
#' # Use all default values
#' model <- ecc(toyml, "RANDOM")
#' pred <- predict(model, toyml)
#'
#' \dontrun{
#' # Use J48 with 100% of instances and only 5 rounds
#' model <- ecc(toyml, 'J48', m = 5, subsample = 1)
#'
#' # Use 75% of attributes
#' model <- ecc(toyml, attr.space = 0.75)
#'
#' # Running in 4 cores and define a specific seed
#' model1 <- ecc(toyml, cores=4, seed=123)
#' }
ecc <- function(mdata, base.method = getOption("utiml.base.method", "SVM"),
                m = 10, subsample = 0.75, attr.space = 0.5, ...,
                cores = getOption("utiml.cores", 1),
                seed = getOption("utiml.seed", NA)) {
  # Validations
  if (class(mdata) != "mldr") {
    stop("First argument must be an mldr object")
  }

  if (m <= 1) {
    stop("The number of iterations (m) must be greater than 1")
  }

  if (subsample < 0.1 || subsample > 1) {
    stop("The subset of training instances must be between 0.1 and 1 inclusive")
  }

  if (attr.space <= 0.1 || attr.space > 1) {
    stop("The attribbute space of training instances must be between 0.1 and 1 inclusive")
  }

  if (cores < 1) {
    stop("Cores must be a positive value")
  }

  # ECC Model class
  eccmodel <- list(rounds = m, call = match.call())
  eccmodel$nrow <- ceiling(mdata$measures$num.instances * subsample)
  eccmodel$ncol <- ceiling(length(mdata$attributesIndexes) * attr.space)


  utiml_preserve_seed()
  if (!anyNA(seed)) {
    set.seed(seed)
  }
  idx <- lapply(seq(m), function(iteration) {
    list(
      rows = sample(mdata$measures$num.instances, eccmodel$nrow),
      cols = sample(mdata$attributesIndexes, eccmodel$ncol),
      chain = sample(rownames(mdata$labels))
    )
  })

  eccmodel$models <- lapply(seq(m), function(iteration) {
    ndata <- create_subset(mdata, idx[[iteration]]$rows, idx[[iteration]]$cols)
    chain <- idx[[iteration]]$chain

    ccmodel <- cc(ndata, base.method, chain, ..., cores = cores, seed = seed)
    ccmodel$attrs <- colnames(ndata$dataset[, ndata$attributesIndexes])
    ccmodel
  })

  utiml_restore_seed()
  class(eccmodel) <- "ECCmodel"
  eccmodel
}

#' Predict Method for Ensemble of Classifier Chains
#'
#' This method predicts values based upon a model trained by \code{\link{ecc}}.
#'
#' @param object Object of class '\code{ECCmodel}'.
#' @param newdata An object containing the new input data. This must be a
#'  matrix, data.frame or a mldr object.
#' @param vote.schema Define the way that ensemble must compute the predictions.
#'  The default valid options are: c("avg", "maj", "max", "min"). If \code{NULL}
#'  then all predictions are returned. (Default: \code{'maj'})
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
#' @seealso \code{\link[=ecc]{Ensemble of Classifier Chains (ECC)}}
#' @export
#'
#' @examples
#' \dontrun{
#' # Predict SVM scores
#' model <- ecc(toyml)
#' pred <- predict(model, toyml)
#'
#' # Predict SVM bipartitions running in 6 cores
#' pred <- predict(model, toyml, probability = FALSE, cores = 6)
#'
#' # Return the classes with the highest score
#' pred <- predict(model, toyml, vote.schema = 'max')
#' }
predict.ECCmodel <- function(object, newdata, vote.schema = "maj",
                             probability = getOption("utiml.use.probs", TRUE),
                             ..., cores = getOption("utiml.cores", 1),
                             seed = getOption("utiml.seed", NA)) {
  # Validations
  if (class(object) != "ECCmodel") {
    stop("First argument must be an ECCmodel object")
  }

  if (cores < 1) {
    stop("Cores must be a positive value")
  }

  utiml_ensemble_check_voteschema(vote.schema)
  utiml_preserve_seed()

  newdata <- utiml_newdata(newdata)
  allpreds <- utiml_lapply(object$models, function(ccmodel) {
    predict.CCmodel(ccmodel, newdata[, ccmodel$attrs], ...)
  }, cores, seed)

  utiml_restore_seed()
  utiml_predict_ensemble(allpreds, vote.schema, probability)
}

#' Print ECC model
#' @param x The ecc model
#' @param ... ignored
#' @export
print.ECCmodel <- function(x, ...) {
    cat("Ensemble of Classifier Chains Model\n\nCall:\n")
    print(x$call)
    cat("\nDetails:")
    cat("\n ", x$rounds, "Iterations")
    cat("\n ", x$nrow, "Instances")
    cat("\n ", x$ncol, "Attributes\n")
}
