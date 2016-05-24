kernParamInit <-
function (kern) {
  
  funcName <- paste(kern$type, "KernParamInit", sep="")
  kern$transforms = list()

  func <- get(funcName, mode="function")
  kern <- func(kern)  

  return (kern)
}
