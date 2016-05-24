
#' Create a MOA classifier
#'
#' Create a MOA classifier
#'
#' @param model character string with a model.
#' E.g. HoeffdingTree, DecisionStump, NaiveBayes, HoeffdingOptionTree, ...
#' The list of known models can be obtained by typing RMOA:::.moaknownmodels. 
#' See the examples and \code{\link{MOAoptions}}.
#' @param control an object of class \code{MOAmodelOptions} as obtained by calling \code{\link{MOAoptions}}
#' @param ... options of parameters passed on to \code{\link{MOAoptions}}, in case \code{control} is left to NULL. 
#' Ignored if \code{control} is supplied
#' @return An object of class \code{MOA_classifier}
#' @seealso \code{\link{MOAoptions}}
#' @export 
#' @examples
#' RMOA:::.moaknownmodels
#' ctrl <- MOAoptions(model = "HoeffdingTree", leafprediction = "MC", 
#'    removePoorAtts = TRUE, binarySplits = TRUE, tieThreshold = 0.20)
#' hdt <- MOA_classifier(model = "HoeffdingTree", control=ctrl)
#' hdt
#' hdt <- MOA_classifier(
#'  model = "HoeffdingTree", 
#'  numericEstimator = "GaussianNumericAttributeClassObserver")
#' hdt
MOA_classifier <- function(model, control=NULL, ...){
  out <- list()
  class(out) <- c(model, "MOA_classifier", "MOA_model")
  out$type <- model
  ## Create the model
  out$moamodel <- .jnew(modelclass(out$type))  
  ## Set MOA options
  if(inherits(control, "MOAmodelOptions")){
    if(control$model != out$type){
      stop("Make control contains options for the correct model")
    }
    out$options <- control
  }else{
    out$options <- MOAoptions(out, ...)   
  }  
  out
}

##' @S3method print MOA_classifier
print.MOA_classifier <- function(x, ...){
  print(x$options)
  try(cat(x$moamodel$toString()), silent=TRUE)
  invisible()
}

#' MOA classification trees
#'
#' MOA classification trees
#'
#' @name MOA_classification_trees
#' @param control an object of class \code{MOAmodelOptions} as obtained by calling \code{\link{MOAoptions}}
#' @param ... options of parameters passed on to \code{\link{MOAoptions}}, in case \code{control} is left to NULL. 
#' Ignored if \code{control} is supplied
#' @return An object of class \code{MOA_classifier} which sets up an untrained MOA model,
#' which can be trained using \code{\link{trainMOA}} 
#' @seealso \code{\link{MOAoptions}}, \code{\link{trainMOA}}
#' @examples
#' ctrl <- MOAoptions(model = "HoeffdingTree", leafprediction = "MC", 
#'    removePoorAtts = TRUE, binarySplits = TRUE, tieThreshold = 0.20)
#' hdt <- HoeffdingTree(control=ctrl)
#' hdt
#' hdt <- HoeffdingTree(numericEstimator = "GaussianNumericAttributeClassObserver")
#' hdt
NULL
#' @export 
#' @rdname MOA_classification_trees
AdaHoeffdingOptionTree <- function(control=NULL, ...) {
  MOA_classifier(model = "AdaHoeffdingOptionTree", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_trees
ASHoeffdingTree <- function(control=NULL, ...) {
  MOA_classifier(model = "ASHoeffdingTree", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_trees
DecisionStump <- function(control=NULL, ...) {
  MOA_classifier(model = "DecisionStump", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_trees
HoeffdingAdaptiveTree <- function(control=NULL, ...) {
  MOA_classifier(model = "HoeffdingAdaptiveTree", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_trees
HoeffdingOptionTree <- function(control=NULL, ...) {
  MOA_classifier(model = "HoeffdingOptionTree", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_trees
HoeffdingTree <- function(control=NULL, ...) {
  MOA_classifier(model = "HoeffdingTree", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_trees
LimAttHoeffdingTree <- function(control=NULL, ...) {
  MOA_classifier(model = "LimAttHoeffdingTree", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_trees
RandomHoeffdingTree <- function(control=NULL, ...) {
  MOA_classifier(model = "RandomHoeffdingTree", control=control, ...)
}


#' MOA bayesian classification
#'
#' MOA bayesian classification
#'
#' @name MOA_classification_bayes
#' @param control an object of class \code{MOAmodelOptions} as obtained by calling \code{\link{MOAoptions}}
#' @param ... options of parameters passed on to \code{\link{MOAoptions}}, in case \code{control} is left to NULL. 
#' Ignored if \code{control} is supplied
#' @return An object of class \code{MOA_classifier} which sets up an untrained MOA model,
#' which can be trained using \code{\link{trainMOA}} 
#' @seealso \code{\link{MOAoptions}}, \code{\link{trainMOA}}
#' @examples
#' ctrl <- MOAoptions(model = "NaiveBayes")
#' mymodel <- NaiveBayes(control=ctrl)
#' mymodel
NULL
#' @export 
#' @rdname MOA_classification_bayes
NaiveBayes <- function(control=NULL, ...) {
  MOA_classifier(model = "NaiveBayes", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_bayes
NaiveBayesMultinomial <- function(control=NULL, ...) {
  MOA_classifier(model = "NaiveBayesMultinomial", control=control, ...)
}

#' MOA active learning classification
#'
#' MOA active learning classification
#'
#' @name MOA_classification_activelearning
#' @param control an object of class \code{MOAmodelOptions} as obtained by calling \code{\link{MOAoptions}}
#' @param ... options of parameters passed on to \code{\link{MOAoptions}}, in case \code{control} is left to NULL. 
#' Ignored if \code{control} is supplied
#' @return An object of class \code{MOA_classifier} which sets up an untrained MOA model,
#' which can be trained using \code{\link{trainMOA}} 
#' @seealso \code{\link{MOAoptions}}, \code{\link{trainMOA}}
#' @examples
#' ctrl <- MOAoptions(model = "ActiveClassifier")
#' mymodel <- ActiveClassifier(control=ctrl)
#' mymodel
NULL
#' @export 
#' @rdname MOA_classification_activelearning
ActiveClassifier <- function(control=NULL, ...) {
  MOA_classifier(model = "ActiveClassifier", control=control, ...)
}


#' MOA classification using ensembles
#'
#' MOA classification using ensembles (bagging/boosting/stacking/other)
#'
#' @name MOA_classification_ensemblelearning
#' @param control an object of class \code{MOAmodelOptions} as obtained by calling \code{\link{MOAoptions}}
#' @param ... options of parameters passed on to \code{\link{MOAoptions}}, in case \code{control} is left to NULL. 
#' Ignored if \code{control} is supplied
#' @return An object of class \code{MOA_classifier} which sets up an untrained MOA model,
#' which can be trained using \code{\link{trainMOA}} 
#' @seealso \code{\link{MOAoptions}}, \code{\link{trainMOA}}
#' @examples
#' ctrl <- MOAoptions(model = "OzaBoostAdwin")
#' mymodel <- OzaBoostAdwin(control=ctrl)
#' mymodel
NULL
#' @export 
#' @rdname MOA_classification_ensemblelearning
AccuracyUpdatedEnsemble <- function(control=NULL, ...) {
  MOA_classifier(model = "AccuracyUpdatedEnsemble", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
AccuracyWeightedEnsemble <- function(control=NULL, ...) {
  MOA_classifier(model = "AccuracyWeightedEnsemble", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
ADACC <- function(control=NULL, ...) {
  MOA_classifier(model = "ADACC", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
DACC <- function(control=NULL, ...) {
  MOA_classifier(model = "DACC", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
LeveragingBag <- function(control=NULL, ...) {
  MOA_classifier(model = "LeveragingBag", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
LimAttClassifier <- function(control=NULL, ...) {
  MOA_classifier(model = "LimAttClassifier", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
OCBoost <- function(control=NULL, ...) {
  MOA_classifier(model = "OCBoost", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
OnlineAccuracyUpdatedEnsemble <- function(control=NULL, ...) {
  MOA_classifier(model = "OnlineAccuracyUpdatedEnsemble", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
OzaBag <- function(control=NULL, ...) {
  MOA_classifier(model = "OzaBag", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
OzaBagAdwin <- function(control=NULL, ...) {
  MOA_classifier(model = "OzaBagAdwin", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
OzaBagASHT <- function(control=NULL, ...) {
  MOA_classifier(model = "OzaBagASHT", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
OzaBoost <- function(control=NULL, ...) {
  MOA_classifier(model = "OzaBoost", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
OzaBoostAdwin <- function(control=NULL, ...) {
  MOA_classifier(model = "OzaBoostAdwin", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
TemporallyAugmentedClassifier <- function(control=NULL, ...) {
  MOA_classifier(model = "TemporallyAugmentedClassifier", control=control, ...)
}
#' @export 
#' @rdname MOA_classification_ensemblelearning
WeightedMajorityAlgorithm <- function(control=NULL, ...) {
  MOA_classifier(model = "WeightedMajorityAlgorithm", control=control, ...)
}




#' Summary statistics of a MOA classifier 
#'
#' Summary statistics of a MOA classifier 
#'
#' @param object an object of class  \code{MOA_classifier}
#' @param ... other arguments, currently not used yet
#' @return the form of the return value depends on the type of MOA model
#' @export 
#' @S3method summary MOA_classifier
#' @examples
#' hdt <- HoeffdingTree(numericEstimator = "GaussianNumericAttributeClassObserver")
#' hdt
#' data(iris)
#' iris <- factorise(iris)
#' irisdatastream <- datastream_dataframe(data=iris)
#' ## Train the model
#' hdtreetrained <- trainMOA(model = hdt, 
#'  Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, 
#'  data = irisdatastream)
#' summary(hdtreetrained$model)
summary.MOA_classifier <- function(object, ...){
  out <- list()
  out$trainingHasStarted <- .jcall(object$moamodel, "Z", "trainingHasStarted")
  out$isRandomizable <- .jcall(object$moamodel, "Z", "isRandomizable")
  out$type <- object$type
  out$options <- object$options$options
  out$fields <- fields(object)[c("attributes","attribute.names","response","responselevels")]
  class(out) <- "summary_MOA_classifier"
  out
}

##' @S3method print summary_MOA_classifier
print.summary_MOA_classifier <- function(x, ...){
  cat(x$type, sep="\n")  
  cat(sprintf("response: %s", x$fields$response), sep="\n")
  cat(sprintf("responselevels: %s", paste(x$fields$responselevels, collapse=", ")), sep="\n")
  cat(sprintf("data features: %s", paste(x$fields$attribute.names, collapse=", ")), sep="\n")
  cat(sprintf("Model has trained: %s", x$trainingHasStarted), sep="\n")
}

