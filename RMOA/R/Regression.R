#' Create a MOA regressor
#'
#' Create a MOA regressor
#'
#' @param model character string with a model.
#' E.g. AMRulesRegressor, FadingTargetMean, FIMTDD, ORTO, Perceptron, RandomRules, SGD, TargetMean, ...
#' The list of known models can be obtained by typing RMOA:::.moaknownmodels. 
#' See the examples and \code{\link{MOAoptions}}.
#' @param control an object of class \code{MOAmodelOptions} as obtained by calling \code{\link{MOAoptions}}
#' @param ... options of parameters passed on to \code{\link{MOAoptions}}, in case \code{control} is left to NULL. 
#' Ignored if \code{control} is supplied
#' @return An object of class \code{MOA_regressor}
#' @seealso \code{\link{MOAoptions}}
#' @export 
#' @examples
#' mymodel <- MOA_regressor(model = "FIMTDD")
#' mymodel
#' data(iris)
#' iris <- factorise(iris)
#' irisdatastream <- datastream_dataframe(data=iris)
#' ## Train the model
#' mytrainedmodel <- trainMOA(model = mymodel, 
#'  Sepal.Length ~ Petal.Length + Species, data = irisdatastream)
#' mytrainedmodel$model
#' 
#' summary(lm(Sepal.Length ~ Petal.Length + Species, data = iris))
#' predict(mytrainedmodel, newdata=iris)
MOA_regressor <- function(model, control=NULL, ...){
  out <- list()
  class(out) <- c(model, "MOA_regressor", "MOA_model")
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

##' @S3method print MOA_regressor
print.MOA_regressor <- function(x, ...){
  print(x$options)
  try(cat(x$moamodel$toString()), silent=TRUE)
  invisible()
}



#' MOA regressors
#'
#' MOA regressors
#'
#' @name MOA_regressors
#' @param control an object of class \code{MOAmodelOptions} as obtained by calling \code{\link{MOAoptions}}
#' @param ... options of parameters passed on to \code{\link{MOAoptions}}, in case \code{control} is left to NULL. 
#' Ignored if \code{control} is supplied
#' @return An object of class \code{MOA_classifier} which sets up an untrained MOA model,
#' which can be trained using \code{\link{trainMOA}} 
#' @seealso \code{\link{MOAoptions}}, \code{\link{trainMOA}}
#' @examples
#' ctrl <- MOAoptions(model = "FIMTDD", DoNotDetectChanges = TRUE, noAnomalyDetection=FALSE,
#'    univariateAnomalyprobabilityThreshold = 0.5, verbosity = 5)
#' mymodel <- FIMTDD(control=ctrl)
#' mymodel
#' mymodel <- FIMTDD(ctrlDoNotDetectChanges = FALSE)
#' mymodel
NULL
#' @export 
#' @rdname MOA_regressors
FIMTDD <- function(control=NULL, ...) {
  MOA_regressor(model = "FIMTDD", control=control, ...)
}
#' @export 
#' @rdname MOA_regressors
AMRulesRegressor <- function(control=NULL, ...) {
  MOA_regressor(model = "AMRulesRegressor", control=control, ...)
}

#' Summary statistics of a MOA regressor 
#'
#' Summary statistics of a MOA regressor 
#'
#' @param object an object of class  \code{MOA_regressor}
#' @param ... other arguments, currently not used yet
#' @return the form of the return value depends on the type of MOA model
#' @export 
#' @S3method summary MOA_regressor
#' @examples
#' ## TODO
summary.MOA_regressor <- function(object, ...){
  out <- list()
  out$trainingHasStarted <- .jcall(object$moamodel, "Z", "trainingHasStarted")
  out$isRandomizable <- .jcall(object$moamodel, "Z", "isRandomizable")
  out$type <- object$type
  out$options <- object$options$options
  out$fields <- fields(object)[c("attributes","attribute.names","response","responselevels")]
  class(out) <- "summary_MOA_regressor"
  out
}

##' @S3method print summary_MOA_regressor
print.summary_MOA_regressor <- function(x, ...){
  cat(x$type, sep="\n")  
  cat(sprintf("response: %s", x$fields$response), sep="\n")
  cat(sprintf("data features: %s", paste(x$fields$attribute.names, collapse=", ")), sep="\n")
  cat(sprintf("Model has trained: %s", x$trainingHasStarted), sep="\n")
}
