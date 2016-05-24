kernDiagCompute <-
function (kern, x) {
  funcName <- paste(kern$type, "KernDiagCompute", sep="")
  func <- get(funcName, mode="function")
  k <- func(kern, x)
  return (k)
}
