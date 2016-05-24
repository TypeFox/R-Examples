print.lordif.MC <-
function(x, ...) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat(paste("  Number of iterations:",x$nr,"\n\n"))
    cat("  Monte Carlo threshold values:\n\n")
    print(x$cutoff)
    invisible(x)
  }
