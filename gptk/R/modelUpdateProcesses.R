modelUpdateProcesses <-
function (model, predt=NULL) {
#   if (is.GPModel(model))
#     return (modelUpdateProcesses(modelStruct(model), predt=predt))

  funcName <- paste(model$type, "UpdateProcesses", sep="")
  func <- get(funcName, mode="function")
  return (func(model, predt=predt))
}
