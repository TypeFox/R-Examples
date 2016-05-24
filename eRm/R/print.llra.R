print.llra <- function(x,...)
  {
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("\nParameters:\n")
    outmat <- rbind(x$etapar,x$se.eta)
    rownames(outmat) <- c("Estimate","Std.Err")
    print(outmat)
  }
