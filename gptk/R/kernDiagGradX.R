kernDiagGradX <-
function (kern, x) {
  funcName <- paste(kern$type, "KernDiagGradX", sep="")
  func <- get(funcName, mode="function")
  k <- func(kern, x)
  return (k)
}
