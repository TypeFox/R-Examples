print.uniscale <- function(x,...)
  {
    cat("\nCall: ")
    print(x$call)
    cat("\nFinal stress value:",x$stress,"\n")
    cat("Number of accepted permutations:",x$npermscale,"\n")
    cat("Number of possible permutations:",x$npermtot,"\n")
    cat("Number of objects:",x$nobj,"\n")
    cat("\n")
  }
