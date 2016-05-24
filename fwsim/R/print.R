print.fwsim <-
function(x, ...) {
  if (!is(x, "fwsim")) stop("x must be a fwsim object")

  if (x$size_model == "stochastic") {
    cat("fwsim (stochastic population size) with ", x$pars$G, " generations with initial population size ", sum(x$pars$N0),
      " and final population size ", sum(x$population[, ncol(x$population)]), ".\n", sep = "")
  } else if (x$size_model == "fixed") {
    cat("fwsim (fixed population size) with ", x$pars$G, " generations with population size ", sum(x$pars$N0), ".\n", sep = "")
  }
  
  return(invisible(NULL))
}

print.mutmodel <-
function(x, ...) {
  if (!is(x, "mutmodel")) stop("x must be a mutmodel object")
  
  if (x$modeltype == 1L) {
    cat("SMM (stepwise mutation model) on ", ncol(x$mutpars), " loci with parameters:\n", sep = "")
  } else if (x$modeltype == 2L) {
    cat("LMM (logistic mutation model) on ", ncol(x$mutpars), " loci with parameters:\n", sep = "")
  } else if (x$modeltype == 3L) {
    cat("EMM (exponential mutation model, unpublished) on ", ncol(x$mutpars), " loci with parameters:\n", sep = "")
  } else {
    stop("Unexpected model type")
  }

  print(x$mutpars)
      
  return(invisible(NULL))
}

