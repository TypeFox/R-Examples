print.genlasso <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  cat("\nOutput:\n")
  cat(paste("Path model with ", length(x$lambda)," total steps.",
            "\n\n", sep=""))
}
