# Create a place for package-wide global variables.
.latentnetEnv <- new.env()

.onLoad <- function(lib, pkg){
  ## Remember where we loaded this instance of latentnet, so that
  ## the snowFT slave functions could do the same.
  .latentnetEnv$path.to.me <- file_path_as_absolute(lib)
  .latentnetEnv$nlog.double.eps <- -log(.Machine[["double.eps"]])
  .latentnetEnv$BIC.warned <- FALSE
}

.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("latentnet",c("statnet"),FALSE)
  if(!is.null(sm)) packageStartupMessage(sm,"NOTE: BIC calculation prior to latentnet 2.7.0 had a bug in the calculation of the effective number of parameters. See help(summary.ergmm) for details.")
}
