modelDisplay <-
function(model, ...) {
#   if (is.GPModel(model))
#     model <- modelStruct(model)

  funcName <- paste(model$type, "Display", sep="")
  if(exists(funcName, mode="function")) {
    func <- get(funcName, mode="function")
    func(model, ...)
  }
}
