print.basta <-
    function(x, ...){
  extraArgs       <- list(...)
  if (length(extraArgs) > 0) {
    if (!is.element('digits', names(extraArgs))){
      digits <- 4
    } else {
      digits <- extraArgs$digits
    }
  } else {
    digits <- 4
  }
  if ("ModelSpecs" %in% names(x)) {
    x$modelSpecs <- x$ModelSpecs
  }
  cat("\nCall:\n")
  cat(paste("Model             \t\t: ", x$modelSpecs[1], "\n", sep = ""))
  cat(paste("Shape             \t\t: ", x$modelSpecs[2], "\n", sep = ""))
  cat(paste("Covars. structure \t\t: ", x$modelSpecs[3], "\n", sep = ""))
  cat(paste("Minimum age       \t\t: ", x$modelSpecs[4], "\n", sep = ""))
  cat(paste("Cat. covars.      \t\t: ", x$modelSpecs[5], "\n", sep = ""))
  cat(paste("Cont. covars.     \t\t: ", x$modelSpecs[6], "\n", 
          collapse = ""))
  
  cat("\nCoefficients:\n")
  print.default(x$coefficients, digits, ...)
  if (x$DIC[1] == "Not calculated"){
    if (x$set['nsim'] == 1) {
      cat("\nConvergence calculations require more than one run.",
          "\nTo estimate potential scale reduction run at least two simulations.\n")
    } else {
      cat("\nWarning: Convergence not reached for some parameters",
          " (i.e. 'PotScaleReduc' values larger than 1.1).",
          "\nThese estimates should not be used for inference.\n")
    }
  } else {
    cat("\nAppropriate convergence reached for all parameters.\n")
    cat("DIC:\n")
    cat(x$DIC["DIC"], "\n")
    if ("Convergence" %in% names(x)) {
      warning("Model fit in versions older than BaSTA 1.5 had a mistake in the",
          "\ncalculation of DIC values. In case you are interested in\n",
          "comparing the fit of different models, please try to run them\n",
          "with BaSTA versions 1.5 or above.",
          " We apologise for the inconvenience.", call. = FALSE)
    }
  } 
}

