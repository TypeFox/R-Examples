kernGradient <-
function (kern, x, ...) {
  funcName <- paste(kern$type, "KernGradient", sep="")
  func <- get(funcName, mode="function")

  g <- func(kern, x, ...)

  factors <- .kernFactors(kern, "gradfact")
  for (i in seq(along=factors))
    g[factors[[i]]$index] <- g[factors[[i]]$index]*factors[[i]]$val

  return (g)
}
