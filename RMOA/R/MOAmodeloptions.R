#' Get and set options for models build with MOA.
#'
#' Get and set options for models build with MOA.
#'
#' @param model character string with a model or an object of class \code{MOA_model}.
#' E.g. HoeffdingTree, DecisionStump, NaiveBayes, HoeffdingOptionTree, ...
#' The list of known models can be obtained by typing RMOA:::.moaknownmodels. See the examples.
#' @param ... other parameters specifying the MOA modelling options of each model. See the examples.
#' @return An object of class \code{MOAmodelOptions}.\cr
#' This is a list with elements:
#'  \enumerate{
#'    \item model: The name of the model
#'    \item moamodelname: The purpose of the model known by MOA (getPurposeString)
#'    \item javaObj: a java reference of MOA options 
#'    \item options: a list with options of the MOA model. Each list element contains the \code{Name} 
#'    of the option, the \code{Purpose} of the option and the current \code{Value}
#'  }
#' See the examples.
#' @export 
#' @examples
#' control <- MOAoptions(model = "HoeffdingTree")
#' control
#' MOAoptions(model = "HoeffdingTree", leafprediction = "MC", 
#'    removePoorAtts = TRUE, binarySplits = TRUE, tieThreshold = 0.20)
#'
#' ## Other models known by RMOA
#' RMOA:::.moaknownmodels
#' 
#' ## Classification Trees
#' MOAoptions(model = "AdaHoeffdingOptionTree")
#' MOAoptions(model = "ASHoeffdingTree")
#' MOAoptions(model = "DecisionStump")
#' MOAoptions(model = "HoeffdingAdaptiveTree")
#' MOAoptions(model = "HoeffdingOptionTree")
#' MOAoptions(model = "HoeffdingTree")
#' MOAoptions(model = "LimAttHoeffdingTree")
#' MOAoptions(model = "RandomHoeffdingTree")
#' ## Classification using Bayes rule
#' MOAoptions(model = "NaiveBayes")
#' MOAoptions(model = "NaiveBayesMultinomial")
#' ## Classification using Active learning
#' MOAoptions(model = "ActiveClassifier")
#' ## Classification using Ensemble learning
#' MOAoptions(model = "AccuracyUpdatedEnsemble")
#' MOAoptions(model = "AccuracyWeightedEnsemble")
#' MOAoptions(model = "ADACC")
#' MOAoptions(model = "DACC")
#' MOAoptions(model = "LeveragingBag")
#' MOAoptions(model = "OCBoost")
#' MOAoptions(model = "OnlineAccuracyUpdatedEnsemble")
#' MOAoptions(model = "OzaBag")
#' MOAoptions(model = "OzaBagAdwin")
#' MOAoptions(model = "OzaBagASHT")
#' MOAoptions(model = "OzaBoost")
#' MOAoptions(model = "OzaBoostAdwin")
#' MOAoptions(model = "TemporallyAugmentedClassifier")
#' MOAoptions(model = "WeightedMajorityAlgorithm")
#' 
#' ## Regressions
#' MOAoptions(model = "AMRulesRegressor")
#' MOAoptions(model = "FadingTargetMean")
#' MOAoptions(model = "FIMTDD")
#' MOAoptions(model = "ORTO")
#' MOAoptions(model = "Perceptron")
#' MOAoptions(model = "SGD")
#' MOAoptions(model = "TargetMean")
MOAoptions <- function(model, ...){
  if(inherits(model, "character")){
    moamodel <- .jnew(modelclass(model))  
  }else{
    moamodel <- model$moamodel
    model <- model$type
  }    
  moaoptions <- moamodel$getOptions()  
  setMOAoptions(moaoptions, ...)
  opts <- list(model = model, 
               moamodelname = .jcall(moamodel, "S", "getPurposeString"), 
               javaObj = moaoptions, 
               options = getMOAoptions(moaoptions))
  class(opts) <- "MOAmodelOptions"
  opts
}


##' @S3method print MOAmodelOptions
print.MOAmodelOptions <- function(x, ...){  
  cat(sprintf("%s modelling options: ", x$model), sep="\n")
  cat(sprintf("MOA model name: %s", x$moamodelname), sep="\n")
  x <- x$options
  for(i in seq_along(x)){
    cat(sprintf("  - %s: %s   (%s)", x[[i]]$Name, x[[i]]$Value, x[[i]]$Purpose), sep="\n")
  }
}


getMOAoptions <- function(x){
  result <- list()
  alloptions <- x$getOptionArray()
  for(i in seq_along(alloptions)){
    optionname <- .jcall(alloptions[[i]], "S", "getName")
    result[[optionname]] <- list()
    result[[optionname]]$Name <- optionname
    result[[optionname]]$Purpose <- .jcall(alloptions[[i]], "S", "getPurpose")    
    result[[optionname]]$Value <- .jcall(alloptions[[i]], "S", "getStateString")
    try(result[[optionname]]$Value <- alloptions[[i]]$getValue(), silent=TRUE)
  }
  result
}

setMOAoptions <- function(x, ...){
  params <- list(...)
  done <- sapply(params, function(x) FALSE)
  
  alloptions <- x$getOptionArray()
  for(i in seq_along(alloptions)){
    optionname <- .jcall(alloptions[[i]], "S", "getName")
    if(optionname %in% names(params)){
      value <- params[[optionname]]
      if(isTRUE(value)){
        value <- tolower(as.character(value))
      }
      .jcall(alloptions[[i]], "V", "setValueViaCLIString", as.character(value))    
      done[[optionname]] <- TRUE
    }    
  }
  if(sum(!done) > 0){
    warning(sprintf("Following MOA options do not exist and are hence not changed: %s", paste(names(done)[done == FALSE], collapse=", ")))
  }
  invisible(done)
}
