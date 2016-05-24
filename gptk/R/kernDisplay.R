kernDisplay <-
function (kern, ...) {

  funcName <- paste(kern$type, "KernDisplay", sep="")
  if(exists(funcName, mode="function")) {
    func <- get(funcName, mode="function")
    return (func(kern, ...))
  }

}
