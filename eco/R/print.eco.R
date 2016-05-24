print.eco <- function(x, digits = max(3, getOption("digits") -3),
                      ...){ 
  cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
  
  if (is.null(x$N))
    N <- rep(1, nrow(x$X))
  else N <- x$N

  W.mean <- cbind(mean(x$W[,1,] %*% (x$X*N/sum(x$X*N))),
                  mean(x$W[,2,] %*% ((1-x$X)*N/sum((1-x$X)*N))))
  colnames(W.mean) <- c("W1", "W2")
  rownames(W.mean) <- "posterior mean"
  
  cat("Aggregate In-sample Estimates:\n\n")
  print.default(format(W.mean, digits = digits), print.gap = 2, quote =
                FALSE)
  cat("\n")
  invisible(x)
}
