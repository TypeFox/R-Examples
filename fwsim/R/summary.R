summary.fwsim <-
function(object, ...) {
  if (!is(object, "fwsim")) stop("object must be a fwsim object")

  if (object$size_model == "stochastic") {
    cat("Fisher-Wright population with stochastic population size.\n", sep = "")
  } else if (object$size_model == "fixed") {
    cat("Fisher-Wright population with fixed population size.\n", sep = "")
  }
  
  cat("Parameters:\n")
  cat("  G        = ", object$pars$G, "\n", sep = "")
  cat("  progress = ", object$pars$progress, "\n", sep = "")
  cat("  trace    = ", object$pars$trace, "\n", sep = "")
  
  if (object$size_model == "stochastic") {
    cat("  alpha    = ", paste(head(object$pars$alpha), collapse = ", "), ", ...\n", sep = "")
  } else if (object$size_model == "fixed") {
    
  }
  
  cat("\n")
  
  cat("H0 and N0:\n")
  print(data.frame(object$pars$H0, N0 = object$pars$N0))
  cat("\n")
  
  #cat("  mutmodel:\n")
  print(object$pars$mutmodel)
  cat("\n")
  
  cat("Number of saved populations = ", sum(unlist(lapply(object$saved_populations, is.matrix))), "\n", sep = "")
  cat("\n")
  
  if (object$size_model == "stochastic") {
    cat("Final population (size = ", sum(object$population[, ncol(object$population)]), "):\n", sep = "")
    print(head(object$population, 5L))
    if (nrow(object$population) > 5L) {
      cat(" (", (nrow(object$population)-5L), " rows hidden.)\n", sep = "")
    } 
    cat("")
    cat("\n")

    cat("Expected population sizes = ", paste(head(object$expected_pop_sizes), collapse = ", "), ", ...\n", sep = "")
    cat("Actual population sizes   = ", paste(head(object$pop_sizes), collapse = ", "), ", ...\n", sep = "")
  } else if (object$size_model == "fixed") {
    cat("Population size = ", sum(object$pars$N0), "\n", sep = "")
  }

  return(invisible(NULL))
}

