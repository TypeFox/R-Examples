modelLogLikelihood <-
function (model) {
#   if (is.GPModel(model))
#     model <- modelStruct(model)

  funcName <- paste(model$type, "LogLikelihood", sep="")
  func <- get(funcName, mode="function")
  return (func(model))
}
