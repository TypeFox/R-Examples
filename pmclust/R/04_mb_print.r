### This files provides functions for output.

mb.print <- function(PARAM, CHECK){
  if(CHECK$algorithm %in% .PMC.CT$algorithm){
    ETA <- PARAM$ETA
    SIGMA <- matrix(do.call("c", PARAM$SIGMA), ncol = PARAM$K)
  }
  MU <- PARAM$MU

  if(.pmclustEnv$COMM.RANK == 0){
    cat("\n")
    cat("Algorithm: ", CHECK$algorithm, "\n", sep = "")
    cat("Convergence: ", CHECK$convergence,
        "  iter: ", CHECK$iter,
        "  abs.err: ", CHECK$abs.err,
        "  rel.err: ", CHECK$rel.err, "\n", sep = "")
    cat("logL: ", PARAM$logL, "\n", sep = "")
    cat("K: ", PARAM$K, "\n", sep = "")
    if(!(CHECK$algorithm %in% c("kmeans", "kmeans.dmat"))){
      cat("\nETA:\n")
      print(ETA)
    }
    cat("\nMU:\n")
    print(MU)
    if(!(CHECK$algorithm %in% c("kmeans", "kmeans.dmat"))){
      cat("\nSIGMA:\n")
      print(SIGMA)
    }
    cat("\n")
  }
} # End of mb.print().
