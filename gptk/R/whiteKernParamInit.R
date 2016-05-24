whiteKernParamInit <-
function (kern) {
  
  kern$variance <- exp(-2)
  kern$nParams <- 1
  kern$paramNames <- c("variance")
  
  kern$transforms <- list(list(index=c(1), type="positive"))

  kern$isStationary <- TRUE

  return (kern)
}
