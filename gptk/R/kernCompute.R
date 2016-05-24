kernCompute <-
function (kern, x, x2) {

  funcName <- paste(kern$type, "KernCompute", sep="")
  func <- get(funcName, mode="function")

  if ( nargs() < 3 ) {
    k <- func(kern, x)
  } else {
    k <- func(kern, x, x2)
  }

  return (k)
}
