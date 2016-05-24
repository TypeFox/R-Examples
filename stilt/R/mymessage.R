mymessage <-
function(optresult) {
  
# Determine if optimization was successful and print a message #!+
  if (optresult$convergence == 0) {
    cat("Optimization SUCCESSFUL!!! Optimization message below:\n\n")
    cat(optresult$message, '\n\n')
    cat("Final parameters\n", names(optresult$par), '\n', optresult$par, '\n\n')
    cat(optresult$iterations, " iterations were performed\n")
    cat("Final likelihood = ", -optresult$objective, '\n\n')
  } else {
    cat("***ERROR*** Optimization FAILED. Optimization message below:\n\n")
    cat(optresult$message, '\n\n')
    cat("If you see a \"Singular Convergence\" message, this indicates that there are probably too\n")
    cat("many parameters. A \"False convergence\" message may mean that the objective\n")
    cat("function gradient was computed incorrectly, or relative tolerance is too\n")
    cat("low, or that either the objective function or its gradient are discontinuous\n")
    cat("near the current iterate of the parameters. For more information on the\n")
    cat("optimization message see PORT documentation: Gay (1990), \"Usage summary for\n")
    cat("selected optimization routines\", available online at\n")
    cat("www.netlib.bell-labs.com/cm/cs/cstr/153.pdf\n\n")
  }

}
