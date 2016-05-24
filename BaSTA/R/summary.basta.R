summary.basta <-
    function(object,...){
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
  if ("ModelSpecs" %in% names(object)) {
    object$modelSpecs <- object$ModelSpecs
  }
  if ("version" %in% names(object)) {
    cat(sprintf("\nOutput from BaSTA version %s\n", object$version))
  }
  cat("\nCall:\n")
  cat(paste("Model             \t\t: ", object$modelSpecs[1], "\n", sep = ""))
  cat(paste("Shape             \t\t: ", object$modelSpecs[2], "\n", sep = ""))
  cat(paste("Covars. structure \t\t: ", object$modelSpecs[3], "\n", sep = ""))
  cat(paste("Minimum age       \t\t: ", object$modelSpecs[4], "\n", sep = ""))
  cat(paste("Cat. covars.      \t\t: ", object$modelSpecs[5], "\n", sep = ""))
  cat(paste("Cont. covars.     \t\t: ", object$modelSpecs[6], "\n", 
          collapse = ""))
  
  cat("\nModel settings:\n")
  print(object$set)
  
  
  cat("\nJumps and priors:\n")
  if ("jumpPriors" %in% names(object)) {
    print(object$jumpPriors)
  } else {
    print(object$JumpP)
  }
  
  cat("\nMean Kullback-Leibler\ndiscrepancy calibration (KLDC):\n")
  if (object$K[1] != "Not calculated") {
    if ("qkl1" %in% names(object$K)) {
      meanKLcalib  <- t((object$K$qkl1 + object$K$qkl2) / 2)
    } else {
      meanKLcalib  <- (object$K$q12 + object$K$q21) / 2
    }
    print.default(meanKLcalib, digits = digits)
  } else {
    if (object$set['nsim'] == 1) {
      cat("KLDC was not calculated due to insufficient number\n",
          " of simulations to estimate convergence.\n")
    } else {
      cat("KLDC was not calculated due to lack of convergence,\n",
          "or because covariates were not included in the model.\n")
    }
  }
  
  
  cat("\nCoefficients:\n")
  print.default(object$coefficients, digits = digits)
  
  cat("\nConvergence:\n")
  if ("Convergence" %in% names(object)) {
    object$convergence <- object$Convergence
  }
  if (object$convergence[1] == "Not calculated") {
    if (object$set['nsim'] == 1) {
      cat("\nConvergence calculations require more than one run.",
          "\nTo estimate potential scale reduction run at least",
          "two simulations.\n")
    } else {
      cat("\nWarning: Convergence not reached for some parameters",
          " (i.e. 'PotScaleReduc' values larger than 1.1).",
          "\nThese estimates should not be used for inference.\n")
    }
  } else {
    if (all(object$convergence[, "Rhat"] < 1.1)) {
      cat("Appropriate convergence reached for all parameters.\n")
    } else {
      cat("\nWarning: Convergence not reached for some parameters",
          " (i.e. 'PotScaleReduc' values larger than 1.1).",
          "\nThese estimates should not be used for inference.\n")
    }
  } 
  cat("\nDIC:\n")
  if (object$DIC[1] != "Not calculated"){
    cat(object$DIC["DIC"],"\n")
    if ("Convergence" %in% names(object)) {
      warning("Model fit in versions older than BaSTA 1.5 had a mistake in the",
          "\ncalculation of DIC values. In case you are interested in\n",
          "comparing the fit of different models, please try to run them\n",
          "with BaSTA versions 1.5 or above.",
          " We apologise for the inconvenience.", call. = FALSE)
    }
  } else {
    if (object$set['nsim'] == 1) {
      cat("DIC was not calculated due to insufficient number",
          "of simulations to estimate convergence.\n")
    } else {
      cat("DIC was not calculated due to lack of convergence.\n")
    }
  }
  ans <- c(list(coefficients = object$coef, DIC = object$DIC,
          KullbackLeibler = object$KullbackLeibler, 
          convergence = object$convergence,
          modelSpecs = object$modelSpecs, settings = object$set))
  return(invisible(ans))
}

