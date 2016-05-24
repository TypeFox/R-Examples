modelObjective <-
function (params, model, ...) {
#   if (any(.packages(all.available=TRUE)=="tigre") && is.GPModel(model))
#     model <- modelStruct(model)

  funcName <- paste(model$type, "Objective", sep="")
  if ( exists(funcName, mode="function") ) {
    func <- get(funcName, mode="function")
    err <- func(params, model, ...)
  } else {
    funcName <- paste(model$type, "ExpandParam", sep="")
    func <- get(funcName, mode="function")
    model <- func(model, params)

    funcName <- paste(model$type, "LogLikelihood", sep="")
    func <- get(funcName, mode="function")
    err <- - func(model, ...)
  }

  return (err)
}
