kernGradX <-
function (kern, x1, x2) {
  funcName <- paste(kern$type, "KernGradX", sep="")
  func <- get(funcName, mode="function")
  k <- func(kern, x1, x2)
  return (k)
}
