print.batwing <-
function(x, ...) {
  cat("Parameters:\n")
  print(x$parameters)
  
  cat("Priors:\n")
  for (i in 1:length(x$priors)) {
    cat("  ", names(x$priors)[i], ": ", paste(x$priors[[i]], collapse = ", "), "\n", sep = "")
  }
   
  cat("Acceptance rate tree:\n")
  print(x$accepted_tree / sum(x$proposals_tree))

  cat("Acceptance rate hyperparameters:\n")
  print(x$accepted_hyperparameters / sum(x$proposals_hyperparameters))

  cat("Last iterations:\n")
  print(tail(x$result))

  invisible(x)
}

