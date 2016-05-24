modelOptimise <-
function (model, options, ...) {
#   if (is.GPModel(model)) {
#     haveGPModel <- TRUE
#     model <- modelStruct(model)
#   }
#   else
#     haveGPModel <- FALSE

  funcName <- paste(model$type, "Optimise", sep="")
  if ( exists(funcName, mode="function") ) {
    func <- get(funcName, mode="function")
    model <- func(model, options, ...)
  } else {
    if ( "optimiser" %in% names(options) ) {
      funcName <- paste(options$optimiser, "optim", sep="")
    } else {
      funcName <- "CGoptim"
    }
    optFunc <- get(funcName, mode="function")

    params <- modelExtractParam(model)
    newParams <- optFunc(params, modelObjective, modelGradient, options, model)
    
    model <- modelExpandParam(model, newParams$xmin)
    model$llscore <- newParams$objective
  }

#   if (haveGPModel)
#     return (new("GPModel", model))
#   else
  return (model)
}
