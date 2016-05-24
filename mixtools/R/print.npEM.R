print.npEM <- function(x, ...) {
  n <- NROW(x$data)
  r <- NCOL(x$data)  
  m <- length(x$lambdahat)
  cat(paste("Observations:", n, "\n"))
  cat(paste("Coordinates per observation:", r, "\n"))
  cat(paste("Mixture components:", m, "\n"))
  if (r>1) {
    B <- max(x$blockid)
    cat(paste("Blocks (of conditionally iid coordinates):", B, "\n\n"))
  }
  dp = match(c("data","posteriors", "lambda", "mu"), names(x), nomatch=0)
  print.default(structure(x[-dp], class=class(x)), ...)
  invisible(x)
}




