rbfKernParamInit <-
function (kern) {
  kern$inverseWidth <- 1
  kern$variance <- 1
  kern$nParams <- 2
  kern$paramNames <- c("inverseWidth", "variance")
  
  kern$isStationary <- TRUE

  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && kern$options$isNormalised)
    kern$isNormalised <- TRUE
  else
    kern$isNormalised <- FALSE

  if ("options" %in% names(kern) && "inverseWidthBounds" %in% names(kern$options)) {
    kern$transforms <- list(list(index=1, type="bounded"),
                            list(index=2, type="positive"))
    kern$transformArgs <- list()
    kern$transformArgs[[1]] <- kern$options$inverseWidthBounds
    kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
  }
  else {
    kern$transforms <- list(list(index=c(1,2), type="positive"))
  }

  return (kern)
}
