### This files provides functions for output.

print.pmclust <- function(x, ...){
  CHECK <- x$check
  PARAM <- x$param
  N.CLASS <- x$n.class

  ETA <- PARAM$ETA
  MU <- PARAM$MU
  SIGMA <- matrix(do.call("c", PARAM$SIGMA), ncol = PARAM$K)

  cat("Algorithm: ", CHECK$algorithm, "\n",
      "Convergence: ", CHECK$convergence,
      "  iter: ", CHECK$iter,
      "  abs.err: ", CHECK$abs.err,
      "  rel.err: ", CHECK$rel.err, "\n", sep = "")
  cat("  N: ", PARAM$N, "  p: ", PARAM$p, "  K: ", PARAM$K,
      "  logL: ", PARAM$logL, "\n", sep = "")
  cat("n.class:", N.CLASS, "\n", sep = " ")
  cat("ETA:\n")
  print(ETA)
  cat("MU: (p by K)\n")
  print(MU)
  # cat("\nSIGMA:\n")
  # print(SIGMA)
  cat("\n")
} # End of print.pmclust().


print.pkmeans <- function(x, ...){
  CHECK <- x$check
  PARAM <- x$param
  N.CLASS <- x$n.class

  MU <- PARAM$MU

  cat("Algorithm: ", CHECK$algorithm, "\n",
      "Convergence: ", CHECK$convergence,
      "  iter: ", CHECK$iter,
      "  abs.err: ", CHECK$abs.err,
      "  rel.err: ", CHECK$rel.err, "\n", sep = "")
  cat("N: ", PARAM$N, "  p: ", PARAM$p, "  K: ", PARAM$K,
      "  logL: ", PARAM$logL, "\n", sep = "")
  cat("n.class:", N.CLASS, "\n", sep = " ")
  cat("MU: (p by K)\n")
  print(MU)
  cat("\n")
} # End of print.pkmeans().
